/*
Copyright (c) 2015 Vladislav Samsonov <vvladxx@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include <iostream>
#include <iomanip>
#include <map>
#include "dsu.h"
#include "triangulation_bucket.h"
#include <SDL2/SDL.h>
#include <GL/gl.h>

static double stopwatch () {
	struct timespec tm;
	clock_gettime(CLOCK_MONOTONIC, &tm);
	static double lasttime = 0;
	double curtime = (double) tm.tv_sec + (double) tm.tv_nsec * 1e-9;
	double diff = curtime - lasttime;
	lasttime = curtime;
	return diff;
}

struct mst_edge {
	double cost;
	size_t i, j;
	
	inline mst_edge (double cost, size_t i, size_t j) : cost(cost), i(i), j(j) {}
	inline bool operator< (const mst_edge & e) const { return cost < e.cost; }
};

typedef Delaunay<double>::point point;

static vector<point> V;
static unsigned vertex_start, vertex_end;
static Delaunay<double> d;
static vector< pair<const point *, const point *> > mst;
static vector<const point *> S;
static vector< pair<const point *, const point *> > M;
static vector< vector<int> > EG;
static vector<size_t> euler_path;
static SDL_Window * window;
static double triangulation_time, mst_time, epm_time, euler_time;

static void reshape (int w, int h) {
	glViewport(0, 0, w, h);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClearDepth(1.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-0.05, 1.05, -0.05, 1.05, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glEnable(GL_BLEND);
	glEnable(GL_ALPHA_TEST);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_MULTISAMPLE);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClear(GL_COLOR_BUFFER_BIT);
}

static void initSDL () {
	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		cerr << "Unable to init SDL, error: " << SDL_GetError() << endl;
		exit(1);
	}
	
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 8);
	
	window = SDL_CreateWindow("Кристофидес", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 800, 600, SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
	
	if (window == nullptr) {
		cerr << "Unable to create SDL window, error: " << SDL_GetError() << endl;
		exit(1);
	}
	
	SDL_GL_CreateContext(window);
	
	reshape(800, 600);
}

static inline void glDrawPoint (float x, float y, double R = 8.0) {
	glPointSize(R);
	glBegin(GL_POINTS);
	glVertex2f(x, y);
	glEnd();
}

static inline void glDrawLine (float x0, float y0, float x1, float y1) {
	glLineWidth(2.5);
	glBegin(GL_LINES);
	glVertex2f(x0, y0);
	glVertex2f(x1, y1);
	glEnd();
}

void input () {
	cin.sync_with_stdio(false);
	cout.sync_with_stdio(false);
	int n;
	cin >> n >> vertex_start >> vertex_end;
	V.resize(n);
	
	for (int i = 0; i < n; ++i) {
		double x, y;
		cin >> x >> y;
		V[i] = point(x, y);
	}
}

static void buildDT () {
	for (auto v : V)
		d.add_point(v);
	d.build();
}

static void buildMST () {
	DSU<int> dsu(V.size());
	vector<mst_edge> edges;
	auto T = d.get_triangles();
	
	for (auto t : T) {
		edges.push_back(mst_edge((* t.a)(* t.b), d.index_by_ptr(t.a), d.index_by_ptr(t.b)));
		edges.push_back(mst_edge((* t.b)(* t.c), d.index_by_ptr(t.b), d.index_by_ptr(t.c)));
		edges.push_back(mst_edge((* t.a)(* t.c), d.index_by_ptr(t.a), d.index_by_ptr(t.c)));
	}
	sort(edges.begin(), edges.end());
	for (auto e : edges) {
		if (dsu.get(e.i) != dsu.get(e.j)) {
			dsu.unite(e.i, e.j);
			mst.push_back(make_pair(d.ptr_by_index(e.i), d.ptr_by_index(e.j)));
		}
	}
}

static void select_wrong_vertex_set () {
	map<const point *, unsigned> degree;
	for (auto e : mst) {
		++degree[e.first];
		++degree[e.second];
	}
	for (auto d : degree) {
		if ((* d.first) == V[vertex_start] || (* d.first) == V[vertex_end]) {
			if (d.second % 2 == 0)
				S.push_back(d.first);
		}
		else {
			if (d.second % 2 == 1)
				S.push_back(d.first);
		}
	}
}

static inline const point * get_nearest_point (const point * p, const vector<const point *> & V) {
	const point * u = V[0];
	double dist = std::numeric_limits<double>::infinity();
	for (auto w : V) {
		if (w == p)
			continue;
		double d = (* w)(* p);
		if (d < dist) {
			u = w;
			dist = d;
		}
	}
	return u;
}

static inline bool admissible (const point & u, const point & v, const point & u_, const point & v_, const double eps) {
	const double dist = u(v);
	return dist <= eps * u(u_) || dist <= eps * v(v_);
}

static void build_perfect_matching (const vector<const point *> & V) {
	const unsigned k = 50 + sqrt(V.size());
	const double eps = (1.0 + (1.0 / k)) * (1.0 + (1.0 / k));
	vector<const point *> W(V);
	for (unsigned i = 0; i < k && W.size(); ++i) {
		const point * v = W[rand() % W.size()];
		const point * u = get_nearest_point(v, W);
		const point * v_ = get_nearest_point(v, V);
		const point * u_ = get_nearest_point(u, V);
		if (admissible(* u, * v, * u_, * v_, eps)) {
			M.push_back(make_pair(u, v));
			auto it = std::find(W.begin(), W.end(), u);
			W.erase(it);
			it = std::find(W.begin(), W.end(), v);
			W.erase(it);
		}
	}
	while (W.size()) {
		const point * v = W[rand() % W.size()];
		const point * u = get_nearest_point(v, W);
		M.push_back(make_pair(u, v));
		auto it = std::find(W.begin(), W.end(), u);
		W.erase(it);
		it = std::find(W.begin(), W.end(), v);
		W.erase(it);
	}
}

static void build_euler_graph () {
	EG.resize(V.size());
	map<point, size_t> point_to_index;
	for (size_t i = 0; i < V.size(); ++i)
		point_to_index[V[i]] = i;
	
	for (auto e : mst) {
		EG[point_to_index[* e.first]].push_back(point_to_index[* e.second]);
		EG[point_to_index[* e.second]].push_back(point_to_index[* e.first]);
	}
	for (auto e : M) {
		EG[point_to_index[* e.first]].push_back(point_to_index[* e.second]);
		EG[point_to_index[* e.second]].push_back(point_to_index[* e.first]);
	}
}

static void find_euler_path (int v) {
	for (size_t i = 0; i < EG[v].size(); ++i) {
		int u = EG[v][i];
		if (u == -1)
			continue;
		EG[v][i] = -1;
		for (size_t j = 0; j < EG[u].size(); ++j)
			if (EG[u][j] == v) {
				EG[u][j] = -1;
				break;
			}
		find_euler_path(u);
	}
	euler_path.push_back(v);
}

static void shortcut_path () {
	vector<bool> visited(V.size());
	reverse(euler_path.begin(), euler_path.end());
	size_t i = 0, j = 0;
	for (; i < euler_path.size(); ++i) {
		if (!visited[euler_path[i]] && euler_path[i] != vertex_end) {
			euler_path[j] = euler_path[i];
			++j;
		}
		visited[euler_path[i]] = true;
	}
	euler_path.resize(j);
	euler_path.push_back(vertex_end);
}

int main () {
	input();
	initSDL();
	
	stopwatch();
	buildDT();
	triangulation_time = stopwatch();
	buildMST();
	mst_time = stopwatch();
	select_wrong_vertex_set();
	stopwatch();
	build_perfect_matching(S);
	epm_time = stopwatch();
	build_euler_graph();
	find_euler_path(vertex_start);
	shortcut_path();
	euler_time = stopwatch();
	
	cout << "Time (triangulation): " << std::setprecision(6) << std::fixed << triangulation_time << 's' << endl;
	cout << "Time (mst): " << mst_time << 's' << endl;
	cout << "Time (perfect matching): " << epm_time << 's' << endl;
	cout << "Time (pathfinding): " << euler_time << 's' << endl;
	
	for (;;) {
		SDL_Event event;
		bool quit = false;
		static unsigned stage = 0;
		const unsigned kStages = 5;
		
		while (SDL_PollEvent(&event)) {
			switch(event.type) {
				case SDL_QUIT:
					quit = true;
					break;
				case SDL_KEYDOWN:
					if (event.key.keysym.sym == SDLK_ESCAPE)
						quit = true;
					else if (event.key.keysym.sym == SDLK_LEFT)
						stage = (kStages + stage - 1) % kStages;
					else if (event.key.keysym.sym == SDLK_RIGHT)
						stage = (kStages + stage + 1) % kStages;
					break;
				case SDL_WINDOWEVENT:
					if (event.window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
						reshape(event.window.data1, event.window.data2);
					break;
			}
		}
	
		if (quit)
			break;
		
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		
		if (stage == 0) {
			auto T = d.get_triangles();
			for (auto t : T) {
				glColor3f(0, 0, 1);
				glDrawLine(t.a->x, t.a->y, t.b->x, t.b->y);
				glDrawLine(t.a->x, t.a->y, t.c->x, t.c->y);
				glDrawLine(t.c->x, t.c->y, t.b->x, t.b->y);
				
			}
		}
		else if (stage == 1) {
			glColor3f(0, 0, 1);
			for (auto e : mst)
				glDrawLine(e.first->x, e.first->y, e.second->x, e.second->y);
		}
		else if (stage == 2) {
			glColor3f(0, 0, 1);
			for (auto e : mst)
				glDrawLine(e.first->x, e.first->y, e.second->x, e.second->y);
			glColor3f(0, 0.8, 0);
			for (auto s : S)
				glDrawPoint(s->x, s->y, 14);
		}
		else if (stage == 3) {
			glColor3f(0, 0, 1);
			for (auto e : mst)
				glDrawLine(e.first->x, e.first->y, e.second->x, e.second->y);
			glColor3f(0, 0.8, 0);
			for (auto e : M)
				glDrawLine(e.first->x, e.first->y, e.second->x, e.second->y);
			glColor3f(0, 0.8, 0);
			for (auto s : S)
				glDrawPoint(s->x, s->y, 14);
		}
		else if (stage == 4) {
			glColor3f(0, 0, 1);
			for (size_t i = 1; i < euler_path.size(); ++i) {
				size_t a = euler_path[i - 1];
				size_t b = euler_path[i];
				glDrawLine(V[a].x, V[a].y, V[b].x, V[b].y);
			}
		}
		
		glColor3f(0, 0, 0);
		for (auto v : V)
			glDrawPoint(v.x, v.y);
		glColor3f(1, 0, 0);
		glDrawPoint(V[vertex_start].x, V[vertex_start].y);
		glDrawPoint(V[vertex_end].x, V[vertex_end].y);
		
		glFlush();
		SDL_GL_SwapWindow(window);
	}
	
	SDL_Quit();
	
	return 0;
}