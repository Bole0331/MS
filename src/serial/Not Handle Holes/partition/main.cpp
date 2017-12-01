#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <sys/time.h>
#include <sys/resource.h>
#include "point.h"

#define next nxt
#define MAX 999999
#define MIN -999999

using namespace std;

typedef pair<pair<int, int>, pair<int, int> > Pair;

vector<point> v;
// for exact correctness 
vector<bool> ok;
vector<int> order;

vector<map<int, int> > next;
priority_queue<pair<double, Pair>, vector<pair<double, Pair> >, greater<pair<double, Pair> > > q;

int times = 0, cnt = 0, parti = 0;

vector<vector<double> > bp;
vector<bool> margin;
vector<int> blockId;

double xmin = MAX, xmax = MIN;
double ymin = MAX, ymax = MIN;
double zmin = MAX, zmax = MIN;

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
	if (x < margin.size() && margin[x]) return;
	order[x] += 1;
	for (map<int, int>::iterator i = next[x].begin(); i != next[x].end(); i++) {
		int a = i->first, b = x;
		if ((a < margin.size() && margin[a])) continue;
		if (a > b) swap(a, b);
		q.push(pair<double, Pair>(value(a, b) + value(b, a), Pair(pair<int, int>(a, order[a]), pair<int,int>(b, order[b]))));
	}
}

bool goodpair(int x, int y) {
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
	int k = v.size(); 
	q.pop();
	v.push_back((v[a] + v[b]) / 2.0); 
	order.push_back(0);
	ok.push_back(true); 
	ok[a] = ok[b] = false; 
	next.push_back(next[a]);
	int l = next[b][a], r = next[a][b];
	next[k].erase(b);
	int t = l;
	while (t != r) {
		next[k][t] = next[b][t];
		t = next[b][t];
	}	
	for (map<int, int>::iterator i = next[k].begin(); i != next[k].end(); i++) {
		int w = i->first;
		for (map<int, int>::iterator j = next[w].begin(); j != next[w].end(); j++) {
			if ((j->first == a || j->first == b) && (j->second != a && j->second != b))
				next[w][k] = j->second;
			if (j->second == a || j->second == b) next[w][j->first] = k;
		}
		next[w].erase(a); 
		next[w].erase(b);
	}
	update(k); 
	for (map<int, int>::iterator i = next[k].begin(); i != next[k].end(); i++) 
		update(i->first);
	return 1;
}

void prepare() {
	for (int i = 0; i < v.size(); i++ ) {
		if (margin[i]) continue;
		for (map<int, int>::iterator j = next[i].begin(); j != next[i].end(); j++ ) 
			if (i < j->first && !margin[j->first]) 
				q.push(pair<double, Pair>(value(i, j->first) + value(j->first, i),
					Pair(pair<int, int>(i, 0),pair<int,int>(j->first, 0))));
	}
}

int getBlockId(point p) {
	int x, y, z;
	for (int i = 1; i <= parti; i++) {
		if (p.x <= bp[0][i]) {
			x = i - 1;
			break;
		}
	}
	for (int i = 1; i <= parti; i++) {
		if (p.y <= bp[1][i]) {
			y = i - 1;
			break;
		}
	}
	for (int i = 1; i <= parti; i++) {
		if (p.z <= bp[2][i]) {
			z = i - 1;
			break;
		}
	}
	return x * parti * parti + y * parti + z;
}

void buildPartition() {
	vector<double> x,y,z;
	double step_x = (xmax - xmin) / parti;
	double step_y = (ymax - ymin) / parti;
	double step_z = (zmax - zmin) / parti;
	x.push_back(xmin); y.push_back(ymin); z.push_back(zmin);
	for (int i = 1; i < parti; i++) {
		x.push_back(xmin+i*step_x);
		y.push_back(ymin+i*step_y);
		z.push_back(zmin+i*step_z);
	}
	x.push_back(xmax); y.push_back(ymax); z.push_back(zmax);
	bp.push_back(x); bp.push_back(y); bp.push_back(z);
	blockId.reserve(v.size());
	for (int i = 0; i < v.size(); i++) {
		blockId.push_back(getBlockId(v[i]));
	}
	for (int i = 0; i < v.size(); i++) {
		if (margin[i]) continue;
		for (map<int, int>::iterator j = next[i].begin(); j != next[i].end(); j++) {
			if (blockId[i] != blockId[j->first]) {
				margin[i] = true;
				margin[j->first] = true;
				break;
			}
		} 
	}
}

void readFile() {
	char flag, head[128];
	while (cin >> flag) {
		if (flag == 'v') {
			point tmp; 
			cin >> tmp.x >> tmp.y >> tmp.z;
			v.push_back(tmp); 
			ok.push_back(true); 
			order.push_back(0);
			margin.push_back(false);
			xmax = max(tmp.x,xmax); xmin = min(tmp.x,xmin);
			ymax = max(tmp.y,ymax); ymin = min(tmp.y,ymin);
			zmax = max(tmp.z,zmax); zmin = min(tmp.z,zmin);
		} else if (flag == 'f') {
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
	buildPartition();
}

int main(int argc, char * argv[]){
	if ( argc != 5 ) {
		printf("usage: ./meshSimplification input output ratio parti\n");
		return 0;
	}
	freopen(argv[1], "r", stdin);
	freopen(argv[2], "w", stdout);
	parti = atoi(argv[4]);
	readFile();
	times = v.size() * (1.0 - atof(argv[3]));
	fprintf(stderr, "Reading file finished\n");
	prepare();
	int avai = 0;
	for (int i = 0; i < v.size(); i++) {
		if (!margin[i]) avai += 1;
	}
	times = min(times,avai - 200);
	int num = 0;
	while(true) {
		num += merge();
		if ( num >= times ) break;
	}
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
