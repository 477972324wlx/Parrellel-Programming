
#include <iostream>
#include <mpi.h>
#include <time.h>
#include <fstream>
#include <queue>
#include <vector>
#include <stack>
#include <algorithm>
#include <unordered_set>
#include <list>

#define NEW_COST_TAG 99
#define INF 0X3F3F3F3F
#define SRC 1


#define debug(x) printf("END HERE?%d\n",x)

using std::vector;
using std::ifstream;
using std::queue;
using std::stack;
using std::list;
typedef std::unordered_set<int> Set;
const int maxn = 120;

int vecSize = 0;
vector<int> vec;


//邻接链表
struct Edge {
	int v;
	int dist;
};


/*
	状态的定义：
	cost : 当前消耗
	pre  : 当前节点
	set	 : 访问过的节点
	vec  : 访问过的节点序列
*/
struct State {
	int cost;
	int pre;
	Set set;
	vector<int> vec;

	void Insert(int x) { vec.push_back(x); set.insert(x); pre = x; }
	bool exist(int x) { return set.count(x); }
	State& operator = (const State& rhs) {
		vec = rhs.vec; pre = rhs.pre;
		set = rhs.set;
		return *this;
	}
	State(vector<int> _vec, int _pre, Set _set, int _cost = 0) :vec(_vec), pre(_pre), set(_set), cost(_cost) {}
	State() { cost = pre = 0; }
};


int nodes, edges;		// 节点，边
vector<Edge> G[maxn];	//邻接链表

bool visit[maxn] = { false };	//是否访问过
int current[maxn] = { false };  //当前弧优化
int backCost[maxn];				//保存 x->src的最小值

State father[maxn];





/*			读入图			*/
void InitGraph(int& nodes, int& edges, vector<Edge> G[]) {

	ifstream fin("in.txt");
	fin >> nodes >> edges;

	int x, y, z;
	for (int i = 0; i < edges; ++i)
	{
		fin >> x >> y >> z;
		G[x].push_back(Edge{ y,z });

		if (y == SRC) {
			backCost[x] = std::min(backCost[x], z);
		}
	}
	fin.close();
}


//初始化搜索树，方法如下
/*
	1. 广度优先搜索，当队列中的状态数目 == 线程数目，进行状态分派
	2. 如果队列中的状态数目 << 线程数目， 直接可以计算
	3. 当队列中的下一次状态迁移时 > 线程数目，最后一个线程分配最后一个状态的父状态
*/
void InitSearchTree(const int& numprocs)
{

	State initState;
	initState.Insert(SRC);

	queue<State> que;
	que.push(initState);

	while (!que.empty()) {
		bool update = false;

		State curState = que.front(); que.pop();
		int pre = curState.pre;

		while (current[pre] < G[pre].size()) {

			Edge nxtEdge = G[pre][current[pre]++];
			int nxtPoint = nxtEdge.v;

			if (curState.exist(nxtPoint)) continue;

			State newState = curState;
			newState.Insert(nxtPoint);
			newState.cost += nxtEdge.dist;
			father[nxtPoint] = curState;	//记录父线程

			que.push(newState);
			update = true;

			if (que.size() >= numprocs) break;
		}
		if (que.size() >= numprocs) break;
		if (!update)break;
	}


	/*
	处理如下：
		1.队列空直接返回
		2.发送序列大小
		3.发送序列
		4.发送cost值
	*/
	for (int i = numprocs - 1; i >= 1; --i) {
		if (que.empty()) return;

		State pState = que.front(); que.pop();
		int siz = int(pState.vec.size());

		MPI_Send(&siz, 1, MPI_INT, i, 666, MPI_COMM_WORLD);

		for (int v : pState.vec)
			MPI_Send(&v, 1, MPI_INT, i, 666, MPI_COMM_WORLD);

		MPI_Send(&pState.cost, 1, MPI_INT, i, 666, MPI_COMM_WORLD);

	}


	/*		给主线程分配工作	*/
	if (!que.empty()) {
		State pState = que.front(); que.pop();
		if (que.empty()) {
			pState = father[pState.pre];
		}
		vecSize = int(pState.vec.size());
		vec = pState.vec;
		vec.push_back(pState.cost);
	}

}



/*
作用：
	初始化每一个线程的搜索初始状态
	记录cost开销
	使用list来记录访问顺序
	*list*的顺序存取效果，头尾插入效果优于vector
*/
void InitSubSearchStage(list<int>& lst, int& src, int& cost) {

	cost = vec[vecSize];
	src = vec[size_t(vecSize) - 1];

	for (int i = 0; i < vecSize; ++i) {
		visit[vec[i]] = true;
		lst.push_back(vec[i]);
	}
}



list<int> minList;			//存储最小的list
int minCostId = -1;			//存储有最小答案的id号
int minCost = 0x3f3f3f3f;   //存储最小值



/*
	单纯本地搜索最小值
	当节点数为nodes时更新
	若当前距离 > minCost 时弹出
*/
void SearchWithoutSending(list<int>& lst, int& src, int& cost, int& nodeCnt) {

	stack<int> stk;
	stack<int> dis;
	stk.push(src);

	while (!stk.empty()) {
		int cur = stk.top();
		if (nodeCnt == nodes) {
			int total_cost = cost + backCost[cur];
			if (total_cost < minCost) {
				minCost = total_cost;
				minList = lst;
			}

			cost -= dis.top();
			dis.pop();

			nodeCnt--;
			current[cur] = 0;
			lst.pop_back();
			visit[cur] = 0;
			stk.pop();
		}
		else {
			bool suc = false;

			while (current[cur] < G[cur].size()) {
				Edge& nextEdge = G[cur][current[cur]++];
				int nextPoint = nextEdge.v;

				if (!visit[nextPoint] && cost + nextEdge.dist < minCost) {
					visit[nextPoint] = 1;
					nodeCnt++;
					cost = cost + nextEdge.dist;
					suc = true;

					stk.push(nextPoint);
					dis.push(nextEdge.dist);
					lst.push_back(nextPoint);
					break;
				}

			}
			if (suc) continue;

			if (!dis.empty()) {
				cost -= dis.top();
				dis.pop();
			}
			current[cur] = 0;
			nodeCnt--;
			visit[cur] = 0;
			lst.pop_back();
			stk.pop();
		}
	}
}




int main(int argc, char* argv[]) {

	memset(backCost, 0x3f, sizeof(backCost));
	backCost[SRC] = 0;


	int myid, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	clock_t start = clock();

	int src, total_cost;
	InitGraph(nodes, edges, G);	//初始化图

	if (myid == 0) {
		//主线程分配
		InitSearchTree(numprocs);
	}
	else {
		int x;
		MPI_Recv(&vecSize, 1, MPI_INT, 0, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//读入大小

		for (int i = 0; i <= vecSize; ++i) {
			MPI_Recv(&x, 1, MPI_INT, 0, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			vec.push_back(x);
			//读入数据 + cost
		}

	};

	int nodeCnt = vecSize;
	list<int> tourList;

	InitSubSearchStage(tourList, src, total_cost);			//初始搜索状态
	SearchWithoutSending(tourList, src, total_cost, nodeCnt);	//	开始搜索

	if (myid != 0) {
		MPI_Send(&minCost, 1, MPI_INT, 0, 666, MPI_COMM_WORLD);
		//子进程发送最小值
	}
	else {
		minCostId = 0;
		int tmpCost = 0;
		for (int i = 1; i < numprocs; ++i) {
			MPI_Recv(&tmpCost, 1, MPI_INT, i, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//主进程存储最小值，存储最小进程id
			if (tmpCost < minCost) {
				minCost = tmpCost;
				minCostId = i;
			}
		}

		for (int i = 1; i < numprocs; ++i) {
			MPI_Send(&minCostId, 1, MPI_INT, i, 666, MPI_COMM_WORLD);
		}
	}

	if (myid != 0) {
		MPI_Recv(&minCostId, 1, MPI_INT, 0, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//发回最小id
	}

	if (minCostId == myid) {
		//如果进程id为最小id，则开始输出
		printf("Used time : %d ms\n", clock() - start);
		if (minCost >= 0x3f3f3f3f) {
			printf("This graph may not be a connected Graph!\n");
		}
		else {

			printf("MinCost = %d\n", minCost);
			for (int v : minList) {
				printf("%d ", v);
			}
			printf("%d\n", SRC);
		}
	}

	MPI_Finalize();

	return 0;
}
