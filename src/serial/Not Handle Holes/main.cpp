#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <sys/time.h>
#include <sys/resource.h>
#include "point.h"

#define next nxt

using namespace std;

typedef pair<pair<int, int>, pair<int, int> > Pair;

vector<point> v;
// for exact correctness 
vector<bool> ok;
vector<int> order;

vector<map<int, int> > next;
priority_queue<pair<double, Pair>, vector<pair<double, Pair> >, greater<pair<double, Pair> > > q;

int times = 0, cnt = 0;

double parallel_time = 0;
long iter_time = 0;
long loop_time = 0;

long vertex = 0;
long triangle = 0;

double error(point &a, point &b, point &c, point &v){
	point p(b[1] * c[2] - b[2] * c[1] + a[2] * c[1] - c[2] * a[1] + a[1] * b[2] - a[2] * b[1],
			b[2] * c[0] - b[0] * c[2] + a[0] * c[2] - a[2] * c[0] + a[2] * b[0] - b[2] * a[0],
			b[0] * c[1] - b[1] * c[0] + a[1] * c[0] - a[0] * c[1] + a[0] * b[1] - a[1] * b[0]);
	p = p.normalize();
	double d = - (p * a);
	d = p * v + d;
	return d * d;
}

double value(int a, int b){
	point mid = (v[a] + v[b]) / 2.0;
	double res = 0;
	for (map<int, int>::iterator j = next[a].begin(); j != next[a].end(); j ++) 
		res += error(v[a], v[j->first], v[j->second], mid);
	return res;
}

void update(int x) {
	order[x] += 1;
	for (map<int, int>::iterator i = next[x].begin(); i != next[x].end(); i++) {
		int a = i->first, b = x;
		if (a > b) swap(a, b);
		q.push(pair<double, Pair>(value(a, b) + value(b, a), Pair(pair<int, int>(a, order[a]), pair<int,int>(b, order[b]))));
	}
}

bool goodpair(int x, int y){
	int cnt = 0;
	for (map<int, int>::iterator i = next[y].begin(); i != next[y].end(); i++)
		if (next[x].find(i->first) != next[x].end()) cnt++;
	return cnt == 2;
}

int merge() {
	int a = q.top().second.first.first, b = q.top().second.second.first;
	if (q.top().second.first.second != order[a] || q.top().second.second.second != order[b]) {
		q.pop(); 
		return 0;
	}
	if (!ok[a] || !ok[b] || (!goodpair(a, b))) {
		q.pop(); 
		return 0;
	}
	fprintf(stderr, "\rFinished: %.2f%%", 100.0*(++cnt) / (times));
	q.pop();
	v[a] = (v[a] + v[b]) / 2.0;
	order[a] = 0;
	ok[b] = false; 
	map<int,int> tmp = next[a];
	int l = next[b][a], r = next[a][b];
	tmp.erase(b);
	int t = l;
	while (t != r) {
		tmp[t] = next[b][t];
		t = next[b][t];
	}	
	for (map<int, int>::iterator i = tmp.begin(); i != tmp.end(); i++) {
		int w = i->first;
		for (map<int, int>::iterator j = next[w].begin(); j != next[w].end(); j++) {
			if ((j->first == a || j->first == b) && (j->second != a && j->second != b))
				next[w][a] = j->second;
			if (j->second == a || j->second == b) j->second = a;
		}
		next[w].erase(b);
	}
	next[a] = tmp;
	update(a); 
	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);
	loop_time += 1;
	for (map<int, int>::iterator i = next[a].begin(); i != next[a].end(); i++) {
		if (i->first >= 0) update(i->first);
		iter_time += 1;
	}
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	parallel_time += span;
	return 1;
}

void readFile() {
	char flag, head[128];
	while (cin >> flag) {
		if (flag == 'v') {
			vertex += 1;
			point tmp; 
			cin >> tmp.x >> tmp.y >> tmp.z;
			v.push_back(tmp); 
			ok.push_back(true); 
			order.push_back(0);
		} else if (flag == 'f') {
			triangle += 1;
			if (next.size() == 0) {
				next.resize(v.size());
				for (int i = 0; i < v.size(); i++ ) 
					next[i].clear();
			}
			int a[5] = {0}; 
			cin >> a[0] >> a[1] >> a[2]; 
			a[3] = a[0], a[4] = a[1];
			for (int i = 0; i < 5; i ++) a[i]--;
			for (int i = 0; i < 3; i ++) next[a[i]][a[i + 1]] = a[i + 2];
		} else cin.getline(head, 128);
	}
	fprintf(stderr, "vertex: %ld triangle: %ld\n", vertex, triangle);
}

void prepare() {
	for (int i = 0; i < v.size(); i++ ) 
		for (map<int, int>::iterator j = next[i].begin(); j != next[i].end(); j++ ) 
			if (i < j->first) 
				q.push(pair<double, Pair>(value(i, j->first) + value(j->first, i),
					Pair(pair<int, int>(i, 0),pair<int,int>(j->first, 0))));
}

int main(int argc, char * argv[]){
	if ( argc != 4 ) {
		printf("usage: ./meshSimplification input ouput ratio\n");
		return 0;
	}
	freopen(argv[1], "r", stdin);
	freopen(argv[2], "w", stdout);
	readFile();
	fprintf(stderr, "Reading file finished\n");
	prepare();
	times = v.size() * (1.0 - atof(argv[3]));
	int num = 0;
	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);
	while(true) {
		num += merge();
		if ( num >= times ) break;
	}
	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	fprintf(stderr, "\ntotal time: %f\n", span);
	fprintf(stderr, "parallel: %f%%\n", parallel_time / span);
	fprintf(stderr, "tasks: %f\n", iter_time * 1.0 / loop_time);
	int n = 0;
	vector<int> id;
	id.resize(v.size());
	for (int i = 0; i < v.size(); i ++){
		if (ok[i]) {
			id[i] = ++ n;
			printf("v %.5lf %.5lf %.5lf\n", v[i].x, v[i].y, v[i].z);
		}
	}
	cerr << "\n";
	for (int i = 0; i < v.size(); i ++) {
		if (!ok[i]) continue;
		for (map<int, int>::iterator j = next[i].begin(); j != next[i].end(); j ++)
			if (j->first > i && j->second > i) 
				printf("f %d %d %d\n", id[i], id[j->first], id[j->second]);
	}
	return 0;
}
