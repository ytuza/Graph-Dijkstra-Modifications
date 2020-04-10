#include <iostream>
#include <vector>
#include <cstddef>
#include <algorithm>
#include <deque>
#include <chrono>
#include <ctime>
#include <thread>    
#include <queue>
#include <fstream>
#include <math.h>  

#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/assume_abstract.hpp>

	
using namespace std;

class CompareDist
{
public:
	bool operator()(pair<int, double> n1, pair<int, double> n2) {
		return n1.second>n2.second;
	}
};

class CompareTreeA
{
public:
	bool operator()(pair<int, double> n1, pair<int, double> n2) {
		return n1.second>n2.second;
	}
};

void pr(vector<int> v) {
	for (int i = 0; i < v.size(); i++)
		cout << v[i] << ",";
	cout << endl;
}

void pr(vector<double> v) {
	for (int i = 0; i < v.size(); i++)
		cout << v[i] << " - ";
	cout << endl;
}

void pr2(vector<pair<int, double>> v) {
	for (int i = 0; i < v.size(); i++)
		cout << "(" << v[i].first << "," << v[i].second << ")-";
	cout << endl;
}

void pr2(vector<pair<int, float>> v) {
	for (int i = 0; i < v.size(); i++)
		cout << v[i].first << ",";
	cout << endl;
}

struct greaterG
{
	template<class T>
	bool operator()(T const &a, T const &b) const { return a.second > b.second; }
};

struct Edge {
	int v, w;
	float peso;
	Edge(int v = -1, int w = -1, float peso = -1) : v(v), w(w), peso(peso) { }
};

class CompareEdge
{
public:
	bool operator()(Edge n1, Edge n2) {
		return n1.peso>n2.peso;
	}
};

void prEdge(vector<Edge> v) {
	float Tpesos = 0;
	for (int i = 0; i < v.size(); i++) {
		//cout << "(" << v[i].v << "," << v[i].w << "),";
		Tpesos = v[i].peso + Tpesos;
		if (i < 15) {
			cout << "peso parcial: " << Tpesos << endl;
		}

	}
	cout << endl;


	cout << "suma de pesos: " << Tpesos << endl;
}


class GraphV
{
	friend class boost::serialization::access;
public:

	GraphV();
	int V() const;
	int E() const;
	vector<pair<double, double>> location;
	double distance(int a, int b) {
		return sqrt(pow(location[a].first - location[b].first,2) +
			pow(location[a].second -location[b].second,2));
	}
	double distance(int a, float x, float y) {
		return sqrt(pow(location[a].first - x, 2) +
			pow(location[a].second - y, 2));
	}

	

	//funciones virtuales puras
	virtual void insert(Edge) = 0;
	virtual void remove(Edge) = 0;
	virtual bool edge(int, int) = 0;
	virtual vector<pair<int, float>> conectados(int u) = 0;
	virtual void edges(vector<Edge> &v) = 0;
	virtual void edgesO(priority_queue<Edge, vector<Edge>, CompareEdge> &v) = 0;
	virtual double cont(int u, int v) = 0;
	
	
	void setDg(bool d) { digraph = d; }
	void setE(int e) { Ecnt = e; }
	void setV(int v) { Vcnt = v; }

	bool getDg(void) { return digraph; }
	vector<int> getCP(void) { return c_points; }
	pair<vector<int>, double> dijkstraG(int s, int fin);
	pair<vector<int>, double> dijkstraX(int s, int fin);
	vector<int> dijkstraW(int s);

	pair<vector<int>, double> TreeA(int s, int fin);
	pair<vector<int>, double> TreeA2(int s, int fin);


	void dfsG(int u, vector<int> &visita);
	void bfsG(int s, vector<bool> &visited);
	void kruskal(vector<Edge> &MST);
	void chargeMR();
	void chargeCP();

	// Ruta mas corta con grilla
	vector<int> caminoG(int a, int fin) {

		vector<int> tem;
		vector<pair<int, float>> u = conectados(a);
		vector<pair<int, float>> v = conectados(fin);

		if (u.empty() || v.empty()) {
			return { tem };
		}

		int x = (location[a].first * 10), y = (location[a].second * 10);
		tem = TreeA2(a, c_points[(y * 10) + x]).first;


		vector<int> prede;
		int j = fin;
		prede.push_back(fin);
		prede.push_back(main_roads[(y * 10) + x][j]);

		
		while (main_roads[(y * 10) + x][j] != c_points[(y * 10) + x]) {
			j = main_roads[(y * 10) + x][j];
			prede.push_back(main_roads[(y * 10) + x][j]);
		}

		int sz = prede.size() - 1;
		for (int x = 0; x < sz; x++, sz--)
		{
			int temp = prede[x];
			prede[x] = prede[sz];
			prede[sz] = temp;
		}
		tem.insert(tem.end(), prede.begin(), prede.end() );

		
		return tem;
	}
	
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		ar & Vcnt;
		ar & Ecnt;
		ar & digraph;
		ar & location;
		ar & c_points;
		ar & main_roads;
	}

	
private:
	int Vcnt, Ecnt;
	bool digraph;
	vector<int> c_points;
	vector<vector<int>> main_roads;

};

BOOST_SERIALIZATION_ASSUME_ABSTRACT(GraphV)

GraphV::GraphV()
{
}

int GraphV::E() const
{
	return Ecnt;
}

int GraphV::V() const
{
	return Vcnt;
}



struct Grupos
{
	int *parent, *rnk;
	int n;

	Grupos(int n)
	{
		this->n = n;
		parent = new int[n + 1];
		rnk = new int[n + 1];

		for (int i = 0; i <= n; i++)
		{
			rnk[i] = 0;
			parent[i] = i;
		}
	}

	int find(int u)
	{

		if (u != parent[u])
			parent[u] = find(parent[u]);
		return parent[u];
	}


	void juntar(int x, int y)
	{
		x = find(x), y = find(y);

		if (rnk[x] > rnk[y])
			parent[y] = x;
		else 
			parent[x] = y;

		if (rnk[x] == rnk[y])
			rnk[y]++;
	}
};


// crea los puntos centrales de cada sector
void GraphV::chargeCP()
{

	vector<double> min(100, 99);
	c_points.resize(100);

	float rango = 0.1;
	for (int i = 0; i < V(); i++) {

		vector<pair<int, float>> a = conectados(i);
		if (a.empty()) continue;

		int x = (location[i].first * 10), y = (location[i].second * 10);
		double tem = distance(i, (x*rango) + 0.05, (y*rango) + 0.05);
		if (tem < min[(y * 10) + x]) {
			min[(y * 10) + x] = tem;
			c_points[(y * 10) + x] = i;
		}
	}
}

// crea los vectores con un dijkstra por cada nodo central
void GraphV::chargeMR()
{
	main_roads.resize(c_points.size());
	for (int i = 0; i < c_points.size(); i++) {
		main_roads[i] = dijkstraW(c_points[i]);
		cout << i << endl;
	}
}


void GraphV::kruskal(vector<Edge> &MST) {
	priority_queue<Edge, vector<Edge>, CompareEdge> cola;
	edgesO(cola);

	cout << V() << " size cola: " << cola.size()<< endl;
	
	float suma = 0;

	Grupos grp(V());

	while (
		//MST.size() < E()-1 &&
		cola.size() > 0) {

		Edge tem = cola.top();
		//cout << " peso: " << tem.peso << endl;
		cola.pop();
		
		int set_u = grp.find(tem.w);
		int set_v = grp.find(tem.v);

		if (set_u != set_v) {
			suma += tem.peso;
			MST.push_back(tem);

			grp.juntar(tem.w, tem.v);
		}
	}
	cout << suma << " pal" << endl;

}
	
// dijkstra sin cola de prioridad
pair<vector<int>, double> GraphV::dijkstraG(int s, int fin) {

	vector<int> pre(V(), -1);
	vector<pair<int, double>> cola;
	vector<double> distancia(V(), 99999);
	vector<bool> visit(V(), false);
	distancia[s] = 0;
	pair<int, double> t(s, 0);

	cola.push_back(t);

	while (!cola.empty()) {
		pair<int, double> tem = cola.back();

		cola.pop_back();

		visit[tem.first] = true;

		vector<pair<int,float>> conec = conectados(tem.first);

		for (int i = 0; i < conec.size(); i++) {
			if (distancia[conec[i].first] > distancia[tem.first] + conec[i].second) {

				distancia[conec[i].first] = distancia[tem.first] + conec[i].second;
				pre[conec[i].first] = tem.first;
				cola.push_back({ conec[i].first, distancia[conec[i].first] });
				sort(cola.begin(), cola.end(), greaterG());
			}
		}
	}

	vector<int> prede;
	int j = fin;
	prede.push_back(fin);
	prede.push_back(pre[j]);

	while (pre[j] != s) {
		j = pre[j];
		prede.push_back(pre[j]);
	}


	int sz = prede.size() - 1;
	for (int x = 0; x < sz; x++, sz--)
	{
		int temp = prede[x];
		prede[x] = prede[sz];
		prede[sz] = temp;
	}

	return { prede,distancia[fin] };
}

// dijkstra con cola de prioridad
pair<vector<int>, double> GraphV::dijkstraX(int s, int fin) {

	vector<int> pre(V(), -1);
	priority_queue<pair<int, double>, vector<pair<int, double>>, CompareDist> cola;

	//vector<pair<int, double>> cola;
	vector<double> distancia(V(), 99999);
	vector<bool> visit(V(), false);
	distancia[s] = 0;
	pair<int, double> t(s, 0);

	cola.push(t);

	while (!cola.empty()) {
		pair<int, double> tem = cola.top();

		cola.pop();

		visit[tem.first] = true;

		vector<pair<int, float>> conec = conectados(tem.first);

		for (int i = 0; i < conec.size(); i++) {
			if (distancia[conec[i].first] > distancia[tem.first] + conec[i].second) {

				distancia[conec[i].first] = distancia[tem.first] + conec[i].second;
				pre[conec[i].first] = tem.first;
				cola.push({ conec[i].first, distancia[conec[i].first] });
				//sort(cola.begin(), cola.end(), greaterG());
			}
		}
	}

	vector<int> prede;
	int j = fin;
	prede.push_back(fin);
	prede.push_back(pre[j]);

	while (pre[j] != s) {
		j = pre[j];
		prede.push_back(pre[j]);
	}


	int sz = prede.size() - 1;
	for (int x = 0; x < sz; x++, sz--)
	{
		int temp = prede[x];
		prede[x] = prede[sz];
		prede[sz] = temp;
	}

	//pr(pre);

	return { prede,distancia[fin] };
}

// SOLO PARA hallar predecesores de todos
vector<int> GraphV::dijkstraW(int s) {

	vector<int> pre(V(), -1);
	priority_queue<pair<int, double>, vector<pair<int, double>>, CompareDist> cola;

	//vector<pair<int, double>> cola;
	vector<double> distancia(V(), 99999);
	vector<bool> visit(V(), false);
	distancia[s] = 0;
	pair<int, double> t(s, 0);

	cola.push(t);

	while (!cola.empty()) {
		pair<int, double> tem = cola.top();

		cola.pop();

		visit[tem.first] = true;

		vector<pair<int, float>> conec = conectados(tem.first);

		for (int i = 0; i < conec.size(); i++) {
			if (distancia[conec[i].first] > distancia[tem.first] + conec[i].second) {

				distancia[conec[i].first] = distancia[tem.first] + conec[i].second;
				pre[conec[i].first] = tem.first;
				cola.push({ conec[i].first, distancia[conec[i].first] });
			}
		}
	}

	return pre;
}


pair<vector<int>, double> GraphV::TreeA(int s, int fin)
{
	vector<int> pre(V(), -1);
	priority_queue<pair<int, double>, vector<pair<int, double>>, CompareDist> cola;

	//vector<pair<int, double>> cola;
	vector<double> distancia(V(), 99999);
	vector<bool> visit(V(), false);
	distancia[s] = 0;
	pair<int, double> t(s, 0);

	cola.push(t);

	while (!cola.empty()) {
		pair<int, double> tem = cola.top();

		cola.pop();

		visit[tem.first] = true;

		vector<pair<int, float>> conec = conectados(tem.first);

		for (int i = 0; i < conec.size(); i++) {
			if (!visit[conec[i].first] && distancia[conec[i].first] > distancia[tem.first] + conec[i].second) {

				distancia[conec[i].first] = distancia[tem.first] + conec[i].second;
				pre[conec[i].first] = tem.first;
				cola.push({ conec[i].first, distancia[conec[i].first] - distance(conec[i].first,fin)});
				//sort(cola.begin(), cola.end(), greaterG());
			}
		}
		if (tem.first == fin) break;
	}

	vector<int> prede;
	int j = fin;
	prede.push_back(fin);
	prede.push_back(pre[j]);

	while (pre[j] != s) {
		j = pre[j];
		prede.push_back(pre[j]);
	}


	int sz = prede.size() - 1;
	for (int x = 0; x < sz; x++, sz--)
	{
		int temp = prede[x];
		prede[x] = prede[sz];
		prede[sz] = temp;
	}

	return { prede,distancia[fin] };

}

pair<vector<int>, double> GraphV::TreeA2(int s, int fin)
{
	vector<int> rpt;
	vector<pair<int, float>> a = conectados(s);
	vector<pair<int, float>> b = conectados(fin);

	if (a.empty() || b.empty()) {
		//cout << "vacio" << endl;
		return {rpt,0};
	}
	vector<pair<int,float>> path;
	vector<bool> visit(V(), false);	
	bool roads = true;
	double edgeD = 0.0;
	path.push_back({ s ,0});
	int next_n;

	while (	path.back().first != fin) {

		pair<int, float> crtNode = path.back();
		float max_d = 99;
		vector<pair<int, float>> conec = conectados(crtNode.first);
		roads = false;
		for (int i = 0; i < conec.size(); i++) {
			if (!visit[conec[i].first] && (distance(conec[i].first,fin)+conec[i].second)<max_d) {
				roads= true;
				edgeD = conec[i].second;
				next_n = conec[i].first;
				max_d = distance(conec[i].first, fin) + conec[i].second;
			}
		}
		if (roads) {
			visit[next_n] = true;
			path.push_back({ next_n , path.back().second + edgeD});
		}
		else {
			path.pop_back();
		}

	}

	//cout << "size path" << path.size() << endl;

	//pr2(path);
	for (int i = 0; i < path.size(); i++) {
		rpt.push_back(path[i].first);
	}

	//cout << "dentro: " << rpt.size()<< " + "<< s << "- "<<fin<<endl;

	return { rpt,path.back().second };

}

void GraphV::dfsG(int u, vector<int> &visita) {

	cout << "(" << u << ") ";
	visita[u] = 1;
	vector<pair<int, float>> conec = conectados(u);

	for (int i = 0; i < conec.size(); i++) {
		if (visita[conec[i].first] == 0) {
			dfsG(conec[i].first, visita);
		}
	}

	visita[u] = 3;
}

void GraphV::bfsG(int s, vector<bool> &visited) {

	vector<int> q;
	visited[s] = true;
	q.push_back(s);
	while (!q.empty())
	{
		s = q.front();
		cout << "(" << s << ") ";
		q.erase(q.begin());

		vector<pair<int, float>> conec = conectados(s);


		for (int i = 0; i<conec.size(); i++)
		{
			if (!visited[conec[i].first])
			{
				visited[conec[i].first] = true;
				q.push_back(conec[i].first);
			}
		}
	}

	pr(q);
}

class GRAPH : public GraphV {
private:
	vector<vector <double> > adj;
	
public:
	GRAPH(int V, bool digraph = false) : adj(V) {
		setDg(digraph);
		setV(V);
		setE(0);

		for (int i = 0; i < V; i++)
			adj[i].assign(V, false);
	}
	bool edge(int, int);
	void insert(Edge);
	void remove(Edge);
	double cont(int u, int v);
	void edges(vector<Edge> &v);
	void edgesO(priority_queue<Edge, vector<Edge>, CompareEdge> &v);


	vector<pair<int, float>> conectados(int u);


	class adjIterator {
		const GRAPH &G;
			int i, v;

	public:
		adjIterator(const GRAPH &G, int v) : G(G), v(v), i(-1) { }

		int beg() { i = -1; return nxt(); }
		int nxt() {
			for (i++; i < G.V(); i++)
				if (G.adj[v][i] != -1) return i;
			return -1;
		}
		bool end() { return i >= G.V(); }
	};
};

vector<pair<int, float>> GRAPH::conectados(int u) {
	vector<pair<int, float>> conec;

	for (int i = 0; i < V(); i++) {
		if (cont(u,i) != 0) {
			conec.push_back({i,cont(u,i)});
		}
	}
	return conec;
}


void GRAPH::edges(vector<Edge> &v) {
	int nE = 0;
	if (v.size() != this->E()) v.resize(this->E());

	for (int i = 0; i < this->V(); i++) {
		typename adjIterator A(*this,i);
		for (int w = A.beg(); !A.end(); w = A.nxt())
			if (this->getDg() || i < w)
				v[nE++] = Edge(i, w, adj[i][w]);
	}
}

void GRAPH::edgesO(priority_queue<Edge, vector<Edge>, CompareEdge> &v) {
	for (int i = 0; i < this->V(); i++) {
		typename adjIterator A(*this, i);
		for (int w = A.beg(); !A.end(); w = A.nxt())
			if (this->getDg() || i < w)
				v.push(Edge(i, w, adj[i][w]));
	}

}

template <class Graph>
void show(const Graph &G) {
	for (int s = 0; s < G.V(); s++) {
		cout.width(2); cout << s << ":";
		typename Graph::adjIterator A(G, s);
		for (int t = A.beg(); !A.end(); t = A.nxt()) {
			cout.width(2); cout << t << " ";
		}
		cout << endl;
	}
}

double GRAPH::cont(int u, int v) {
	return adj[u][v];
}

void GRAPH::insert(Edge e)
{
	adj[e.v][e.w] = e.peso;
	setE(E() + 1);
}

void GRAPH::remove(Edge e)
{
	adj[e.v][e.w] = 0;
	setE(E() - 1);
}

bool GRAPH::edge(int a, int b)
{
	return adj[a][b];
}

class GraphList: public GraphV
{
	friend class boost::serialization::access;
public:
	GraphList() {}
	GraphList(int V, bool digraph = false) : adj(V){
		setDg(digraph);
		setV(V);
		setE(0);
	}
	~GraphList();

	void insert(Edge e) {
		/*
		for (int i = 0; i < adj[e.w].size(); i++) {
			if (adj[e.w][i].first == e.v) return;
		}
		*/

		adj[e.v].push_back({e.w,e.peso});
		setE(E() + 1);
	}
	void remove(Edge e) {

		vector<pair<int, float>>::iterator it = adj[e.v].begin();


		while ((*it).first != e.w && it != adj[e.v].end()) {
			it++;
		}

		if ((*it).first == e.w) {
			adj[e.v].erase(it);
		}
		
		setE(E() - 1);
	}
	bool edge(int a, int b) {
		for (int i = 0; i < adj[a].size(); i++)
			if (adj[a][i].first == b)
				return true;
		return false;
	}
	
	void edges(vector<Edge> &v) {
		int E = 0;
		if (v.size() != this->E()) v.resize(this->E());
		for (int i = 0; i < this->V(); i++) {
			vector<pair<int, float>> conec = conectados(i);

			for (int j = 0; j < conec.size(); j++)
				if (this->getDg() || i < conec[j].first)
					v[E++] = Edge(i, conec[j].first, conec[j].second);
		}
	}

	void edgesO(priority_queue<Edge, vector<Edge>, CompareEdge> &v) {
		for (int i = 0; i < this->V(); i++) {
			vector<pair<int, float>> conec = conectados(i);
			for (int j = 0; j < conec.size(); j++)
				if (this->getDg() || i < conec[j].first)
					v.push(Edge(i, conec[j].first, conec[j].second));
		}

	}

	vector<pair<int, float>> conectados(int u) {
		return adj[u];
	}

	double cont(int u, int v) {
		vector<pair<int, float>> p = conectados(u);
		for (int j = 0; j < p.size(); j++) {
			if (p[j].first == v) return p[j].second;
		}
		return 0;
	}

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version)
	{
		// save/load base class information
		ar & boost::serialization::base_object<GraphV>(*this);
		ar & adj;
	}

private:
	vector<vector<pair<int, float>>> adj;

};



void save_gl(const GraphList &s) {
	ofstream file{ "10MgPre.bin", ios::binary | ios::trunc };
	boost::archive::binary_oarchive oa{ file };

	oa << s;

	file.flush();
	file.close();
}

void restore_gl(GraphList &s)
{
	ifstream file{ "10MgPre.bin", ios::binary };
	boost::archive::binary_iarchive ia{ file };

	ia >> s;
	file.close();
}

GraphList:: ~GraphList()
{
}




//Dijkstra sin cola de prioridad y solo con matriz de adyacencia
pair<vector<int>, double> dijkstra(GRAPH &G, int s,int fin) {

	cout << "nodos " << G.V() << endl;
	vector<int> pre(G.V(),-1);
	vector<pair<int, double>> cola;
	vector<double> distancia(G.V(), 99999);
	vector<bool> visit(G.V(), false);
	distancia[s]=0;
	pair<int, double> t(s, 0);

	cola.push_back(t);

	int b;
	while (!cola.empty()) {
		//cout << "cola T: " << cola.size()<< endl;
		pair<int, double> tem = cola.back();

		cola.pop_back();

		visit[tem.first] = true;


		for (int i = 0; i < G.V(); i++) {
			if (G.cont(tem.first, i) != 0 && 
				distancia[i] > distancia[tem.first] + G.cont(tem.first,i) ) {

				//cout << "i = " << i << endl;
				distancia[i] = distancia[tem.first] + G.cont(tem.first, i);
				pre[i] = tem.first;
				cola.push_back({i, distancia[i]});
				sort(cola.begin(), cola.end(), greaterG());
			}
		}	

		//pr2(cola);
	}

	vector<int> prede;
	int j = fin;
	prede.push_back(fin);
	prede.push_back(pre[j]);
	
	//cout << "termino el ciclo" << endl;

	while (pre[j] != s) {
		j = pre[j];
		prede.push_back(pre[j]);
		//cout << "tamanio: "<< prede.size() << "  pre: "<< pre[j] <<endl;
	}


	int sz = prede.size()-1;
	for (int x = 0; x < sz; x++, sz--)
	{
		int temp = prede[x];
		prede[x] = prede[sz];
		prede[sz] = temp;
	}

	//cout << "termino el djk" << endl;

	return { prede,distancia[fin] };
}


void dT(GraphList &a,int z, int b) {
	a.dijkstraX(z,b);
}

void linealAT(GraphList &gr, int rep) {
	srand(time(NULL));
	using namespace std::chrono;

	double promedio = 0;
	for (int i = 0; i < rep; i++) {
		int a = rand() % 1000000, b = rand() % 1000000;
		//cout << a << " - " << b << endl;
		//pr(gr.caminoG(a, b));
		gr.caminoG(a, b); // implementada en la clase virtual 
	}
}

void hilosAT(GraphList &a, int total) {
	unsigned int n = std::thread::hardware_concurrency();
	int s = total / n;
	//cout << "lanzando " << n << " threads con " << s << " iteraciones cada uno" << endl;
	vector<thread> vt;

	for (int i = 0; i < n - 1; i++) {
		vt.push_back(thread(linealAT,a,s));
	}
	for (auto &t : vt) {
		t.join();
	}
}

void pruebaLAB(GraphList &a) {
	using namespace std::chrono;

	for(int i = 10 ; i != 10000 ;i=i*10){
		cout << "Nro deconsultas: " << i << endl;
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		linealAT(a, i);
		high_resolution_clock::time_point t2 = high_resolution_clock::now();

		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

		std::cout << "Lineal duro " << time_span.count() << " segundos." << endl;

		t1 = high_resolution_clock::now();
		hilosAT(a, i);
		t2 = high_resolution_clock::now();

		time_span = duration_cast<duration<double>>(t2 - t1);

		std::cout << "Con threads duro " << time_span.count() << " segundos." << endl;
	}
}

int main(){

	/*
	cout << "comenzo" << endl;
	int mV ,mE;
	cin >> mV >> mE;
	GraphList grafo(mV,mE);

	for (int i = 0; i < mV; i++) {
		double a, b;
		cin >> a >> b;
		grafo.location.push_back({ a,b });

	}

	for (int i = 0; i < mE; i++) {
		
		int a, b;
		double w;
		cin >> a >> b >> w;
		Edge e(a, b, w);
		grafo.insert(e);

	}

	
	
	cout << "nodos :" << mV << " location: " << grafo.location.size() << endl;
	
	cout << "termino de llenar" << endl;
	

	grafo.chargeCP();
	//termino puntos
	grafo.chargeMR();
	cout << "pre cal" << endl;

	save_gl(grafo);
	
	*/
		
	GraphList gr;
	restore_gl(gr);
	pruebaLAB(gr);

	/*
	using namespace std::chrono;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	restore_gl(gr);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

	std::cout << "el tiempo de carga es " << time_span.count() << " seconds." << endl;
	cout << "vectore" << gr.V() << endl;
	



	int p = 5;
	t1 = high_resolution_clock::now();
	//hilosAT(gr,10000);
	linealAT(gr, p);
	//gr.TreeA2(10, 100);
	//gr.caminoG(10, 100);
	t2 = high_resolution_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	*/
	/*
	
	std::cout << "se ejecuto en " << time_span.count() << " seconds." << endl;
	t1 = high_resolution_clock::now();
	hilosAT(gr,p);
	//linealAT(gr, 10);
	//gr.TreeA2(10, 100);
	//gr.caminoG(10, 100);
	t2 = high_resolution_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);

	std::cout << "se ejecuto en " << time_span.count() << " seconds." << endl;
	*/
	
	return 0;
}