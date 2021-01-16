import math
import random
import sys
import copy
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

class HybridAntColony:
    
    # 各类参数初始化，接受z矩阵
    # z矩阵每条数据格式应为
    # [0,[12.50,37.60],600,[9.5,10.5],1/6]
    def __init__(self, z):

        super().__init__()
        self.F = 200                    #冷藏车固定成本
        self.F_EMPTY_RATIO = 0.18       #空车燃料消耗率
        self.F_FULL_RATIO = 0.41        #满载燃料消耗率
        self.GREEN_RATIO = 5.49         #绿色成本系数，即柴油价加税
        self.Q_MAX = 9000               #车辆最大载重
        self.A = 5                      #运输消耗系数
        self.B = 12                     #装卸消耗系数
        self.P = 12                     #单位价格
        self.SPEED = 40                 #行驶速度
        self.DEPART_TIME = 6.5            #出发时间
        self.NOW_TIME = 6.5               #当前时间
        self.TRANS_REDUCT_RATIO = 0.0025#运输衰减率
        self.LOAD_REDUCT_RATIO = 0.005  #装卸衰减率
        self.FORWARD_PUNISH_RATIO = 20  #提前惩罚系数
        self.LATE_PUNISHI_RATIO = 20   #迟到惩罚系数
        self.alpha = 1                  #信息素因子
        self.beta = 3                   #启发素因子
        self.sita = 0.8                 #信息素挥发因子
        self.q_0 = 0.6                  #参数q0
        self.Q = 100                    #信息素总量
        self.ant = 15                   #蚂蚁数量
        self.lambda_phe = 1.5           #lambda取1.8
        self.eta = 0.99                 #eta取0.99
        self.times = 100                #迭代次数为100
        self.N_MAX = 5
        self.C_FRESH = 0.4                #新引入参数
        self.gamma = 5                  #朱申俊引入的beta

        self.best_solution = Solution() #全局最优方案
        self.now_best = Solution()      #当前最优方案
        self.best_solution.totalcost = sys.maxsize
        self.r = 0                      #当前节点
        self.customers = []             #客户
        self.graph = []                 #邻接矩阵
        self.untreated = []             #未处理客户节点
        self.solution = Solution()      #方案
        self.pheromone = []             #信息素矩阵
        self.herustic = []              #启发素矩阵
    
        # 将结点信息转化为类
        for i in range(0,len(z)):
            self.customers.append(Customer())
            self.customers[i].id = z[i][0]
            self.customers[i].location_x = z[i][1][0]
            self.customers[i].location_y = z[i][1][1]
            self.customers[i].q_i = z[i][2]
            self.customers[i].e_i = z[i][3][0]
            self.customers[i].l_i = z[i][3][1]
            self.customers[i].t_load = z[i][4]

        # 邻接矩阵记录各节点之间的距离，初始化信息素和启发式因子矩阵
        for i in range(0,len(self.customers)):
            self.graph.append([])
            self.pheromone.append([])
            self.herustic.append([])
            for j in range(0,len(self.customers)):
                self.pheromone[i].append([])
                self.herustic[i].append([])
                self.graph[i].append(math.sqrt(pow((self.customers[j].location_y - self.customers[i].location_y),2) + pow((self.customers[j].location_x - self.customers[i].location_x),2)))

    # 信息素第一次初始化
    def init(self):

        num = len(self.customers) * (len(self.customers)-1) / 2
        phe_0 = self.Q / (2 * num)

        for i in range(0,len(self.customers)-1):
            for j in range(i+1,len(self.customers)):
                self.pheromone[i][j] = self.pheromone[j][i] = phe_0
                # 按照朱申俊方法设计启发式因子
                self.herustic[i][j] = self.herustic[j][i] = self.customers[j].q_i * (1 / self.graph[i][j]) * (1 / (self.customers[j].t_load * 60))     # * 60
    
    # 初始化未服务客户节点和当前节点
    def reset(self):

        # 未服务客户节点不包括配送中心，长度比customers小1
        for i in range(0,len(self.customers)-1):
            self.untreated.append(i+1)
        
        self.solution = Solution()
        self.r = 0
        self.NOW_TIME = 6.5

    def a_star(self):

        # 初始化一条路线
        route = Route()
        # 路线的customers字段添加0
        route.customers.append(0)
        
        # 可能选为下一个节点的中间变量
        min_untreated = 0

        while len(self.untreated) != 0:
            # 初始化成本变量
            C = sys.maxsize
            for i in self.untreated:
                self.customers[i].t_i = self.NOW_TIME + self.graph[self.r][i] / self.SPEED
                c_2 = (self.F_EMPTY_RATIO + ((self.F_FULL_RATIO - self.F_EMPTY_RATIO) * route.load / self.Q_MAX)) * self.graph[self.r][i] * self.GREEN_RATIO
                c_3 = self.A * self.graph[self.r][i] / self.SPEED + self.B * self.customers[i].t_load
                c_41 = self.P * (self.customers[i].q_i * (1 - math.exp(-self.TRANS_REDUCT_RATIO * (self.customers[i].t_i - 6.5))))
                c_42 = self.P * route.load * (1 - math.exp(-self.LOAD_REDUCT_RATIO * self.customers[i].t_load))
                c_4 = c_41 + c_42
                c_5 = 0.0
                if self.customers[i].e_i - self.customers[i].t_i > 0:
                    c_5 = self.FORWARD_PUNISH_RATIO * (self.customers[i].e_i - self.customers[i].t_i)
                elif self.customers[i].t_i - self.customers[i].l_i > 0:
                    c_5 = self.LATE_PUNISHI_RATIO * (self.customers[i].t_i - self.customers[i].l_i)
                c_piece = c_2 + c_3 + c_4 + c_5
                if c_piece < C:
                    C = c_piece
                    min_untreated = i
            load = route.load + self.customers[min_untreated].q_i
            if load > self.Q_MAX:
                route.customers.append(0)
                self.NOW_TIME = 6.5
                self.best_solution.routes.append(route)
                route = Route()
                route.customers.append(0)
                self.r = 0
            route.load += self.customers[min_untreated].q_i
            route.customers.append(min_untreated)
            self.NOW_TIME += (self.customers[min_untreated].t_load + (self.graph[self.r][min_untreated] / self.SPEED))
            self.r = min_untreated
            for i in self.untreated:
                if i == min_untreated:
                    self.untreated.remove(i)
                    break
        route.customers.append(0)
        self.best_solution.routes.append(route)
        self.best_solution.c1 = len(self.best_solution.routes) * self.F
        for i in self.best_solution.routes:
            for j in range(0,len(i.customers)-1):
                self.best_solution.c2 += (self.F_EMPTY_RATIO + (self.F_FULL_RATIO - self.F_EMPTY_RATIO) * i.load / self.Q_MAX) * self.graph[i.customers[j]][i.customers[j+1]] * self.GREEN_RATIO
                self.best_solution.c3 += self.A * (self.graph[i.customers[j]][i.customers[j+1]] / self.SPEED) + self.B * self.customers[i.customers[j+1]].t_load
                #新加的
                self.best_solution.c3 += self.C_FRESH * (i.load / 100) * ((self.graph[i.customers[j]][i.customers[j+1]] / self.SPEED) + self.customers[i.customers[j+1]].t_load)
                #有改动
                self.best_solution.c4 +=self.P * ((self.customers[i.customers[j+1]].q_i * (1 - math.exp(-(self.TRANS_REDUCT_RATIO / (self.gamma * self.C_FRESH)) * (self.customers[i.customers[j+1]].t_i - 6.5)))) + (i.load - self.customers[i.customers[j+1]].q_i) * (1 - math.exp(-(self.LOAD_REDUCT_RATIO / (self.gamma * self.C_FRESH)) * self.customers[i.customers[j+1]].t_load)))
                if (self.customers[i.customers[j+1]].e_i - self.customers[i.customers[j+1]].t_i) > 0:
                    self.best_solution.c5 += self.FORWARD_PUNISH_RATIO * (self.customers[i.customers[j+1]].e_i - self.customers[i.customers[j+1]].t_i)
                elif (self.customers[i.customers[j+1]].t_i - self.customers[i.customers[j+1]].l_i) > 0:
                    self.best_solution.c5 += self.LATE_PUNISHI_RATIO * (self.customers[i.customers[j+1]].t_i - self.customers[i.customers[j+1]].l_i)
                i.load -= self.customers[i.customers[j+1]].q_i
        self.best_solution.totalcost = self.best_solution.c1 + self.best_solution.c2 + self.best_solution.c3 + self.best_solution.c4 + self.best_solution.c5

    # 信息素第二次初始化
    def init_phe(self):

        for i in self.best_solution.routes:
            for j in range(0,len(i.customers)-1):
                self.pheromone[i.customers[j]][i.customers[j+1]] *= self.lambda_phe
                self.pheromone[i.customers[j+1]][i.customers[j]] = self.pheromone[i.customers[j]][i.customers[j+1]]

    # 创建方案
    def construct_solution(self):

        route = Route()
        route.customers.append(0)

        while len(self.untreated) != 0:
            next_node = self.select_next(route)
            if next_node == 0:
                route.customers.append(0)
                # 路线总时间
                route.time += self.graph[self.r][0] / self.SPEED
                self.NOW_TIME += self.graph[self.r][0] / self.SPEED
                # 路线总长度
                route.distance += self.graph[self.r][0]
                # 当前方案添加该行车路线
                self.solution.routes.append(route)
                # 新建一条行车路线，方法同前面
                route = Route()
                route.customers.append(0)
                self.r = 0
                self.NOW_TIME = 6.5
            else:
                # 路线添加客户
                route.customers.append(next_node)
                self.customers[next_node].t_i = self.NOW_TIME + self.graph[self.r][next_node] / self.SPEED
                self.solution.c3 += self.A * (self.graph[self.r][next_node] / self.SPEED) + self.B * self.customers[next_node].t_load
                # 路线总载重添加客户节点的需求
                route.load += self.customers[next_node].q_i
                # 路线总时间添加
                route.time += self.graph[self.r][next_node] / self.SPEED + self.customers[next_node].t_load
                self.NOW_TIME += self.graph[self.r][next_node] / self.SPEED + self.customers[next_node].t_load
                # 路线总长度
                route.distance += self.graph[self.r][next_node]
                # 更改当前节点
                self.r = next_node
                # 将未服务列表移除当前节点
                for i in self.untreated:
                    if i == next_node:
                        self.untreated.remove(i)
                        break

        # 最后还要添加0节点
        route.customers.append(0) 
        # 路线总时间       
        # route.time = 
        # 路线总长度
        route.distance += self.graph[self.r][0]
        # 当前方案添加该行车路线
        self.solution.routes.append(route)
        # 当前方案总消费
        self.solution.c1 = len(self.solution.routes) * self.F
        for i in self.solution.routes:
            for j in range(0,len(i.customers)-1):
                self.solution.c2 += (self.F_EMPTY_RATIO + (self.F_FULL_RATIO - self.F_EMPTY_RATIO) * i.load / self.Q_MAX) * self.graph[i.customers[j]][i.customers[j+1]] * self.GREEN_RATIO
                #新加的
                self.solution.c3 += self.C_FRESH * (i.load / 100) * ((self.graph[i.customers[j]][i.customers[j+1]] / self.SPEED) + self.customers[i.customers[j+1]].t_load)
                #有改动
                self.solution.c4 += self.P * ((self.customers[i.customers[j+1]].q_i * (1 - math.exp(-(self.TRANS_REDUCT_RATIO / (self.gamma * self.C_FRESH)) * (self.customers[i.customers[j+1]].t_i - 6.5)))) + (i.load - self.customers[i.customers[j+1]].q_i) * (1 - math.exp(-(self.LOAD_REDUCT_RATIO / (self.gamma * self.C_FRESH)) * self.customers[i.customers[j+1]].t_load)))
                if (self.customers[i.customers[j+1]].e_i - self.customers[i.customers[j+1]].t_i) > 0:
                    self.solution.c5 += self.FORWARD_PUNISH_RATIO * (self.customers[i.customers[j+1]].e_i - self.customers[i.customers[j+1]].t_i)
                elif (self.customers[i.customers[j+1]].t_i - self.customers[i.customers[j+1]].l_i) > 0:
                    self.solution.c5 += self.LATE_PUNISHI_RATIO * (self.customers[i.customers[j+1]].t_i - self.customers[i.customers[j+1]].l_i)
                i.load -= self.customers[i.customers[j+1]].q_i
        self.solution.totalcost = self.solution.c1 + self.solution.c2 + self.solution.c3 + self.solution.c4 + self.solution.c5

    # 选择下一个节点，接收参数route:当前路线
    def select_next(self, route):

        # 如果未服务列表为空，返回0点
        if len(self.untreated) == 0:
            return 0
        
        # 初始化下一个节点和随机的q1
        next_node = 0
        random_num = random.random()
        # 如果q0 >= q1
        if self.q_0 >= random_num:
            # 初始化最大信息素和启发因子的积（带参数α和β）
            max_phe = 0
            max_untreated_num = 0
            # 遍历未服务节点
            for i in range(0,len(self.untreated)):
                # 找最大信息素和启发因子的积（带参数α和β）
                if math.pow(self.pheromone[self.r][self.untreated[i]], self.alpha) * math.pow(self.herustic[self.r][self.untreated[i]], self.beta) > max_phe:
                    max_phe = math.pow(self.pheromone[self.r][self.untreated[i]], self.alpha) * math.pow(self.herustic[self.r][self.untreated[i]], self.beta)
                    max_untreated_num = i
            # 找到最大积对应的点设为下一个节点
            next_node = self.untreated[max_untreated_num]
            # 判断是否超载，超载则返回0
            load = route.load + self.customers[next_node].q_i
            if load > self.Q_MAX:
                next_node = 0
        # q0 < q1
        else:
            # 累加概率初始化
            sum_p = 0
            # 概率区间列表
            p_i = []
            # 信息素和启发因子的积的总和
            sum_phe_her = 0
            for i in range(0,len(self.untreated)):
                sum_phe_her += math.pow(self.pheromone[self.r][self.untreated[i]], self.alpha) * math.pow(self.herustic[self.r][self.untreated[i]], self.beta)
            for i in range(0,len(self.untreated)):
                p = math.pow(self.pheromone[self.r][self.untreated[i]], self.alpha) * math.pow(self.herustic[self.r][self.untreated[i]], self.beta) / sum_phe_her
                p_i.append([sum_p, sum_p+p])
                sum_p += p
            # 轮盘赌随机数
            random_num_2 = random.random()
            # 判断随机数落在哪个概率区间，落在哪里就将对应的未服务节点设为下一个节点
            for i in range(0,len(p_i)):
                if random_num_2 >= p_i[i][0] and random_num_2 < p_i[i][1]:
                    next_node = self.untreated[i]
        
        # 再次验证是否超载，如超载，将返回0点
        load = route.load + self.customers[next_node].q_i
        if load > self.Q_MAX:
            next_node = 0

        return next_node

    # 更新信息素,要改一些变量
    def update_pheromone(self, solutions, delta_max, delta_min):

        # 初始化信息素增加值
        delta = 0

        ant_solution_node = []
        ant_solution_path = []
        num = 0
        self.now_best.totalcost = sys.maxsize
        for i in solutions:
            ant_solution_node.append([])
            if i.totalcost < self.now_best.totalcost:
                self.now_best = i
            for j in i.routes:
                for k in j.customers[:-1]:
                    ant_solution_node[num].append(k)
            ant_solution_node[num].append(0)
            num += 1

        num = 0
        for i in ant_solution_node:
            ant_solution_path.append([])
            for j in range(0,len(i)-1):
                ant_solution_path[num].append([i[j],i[j+1]])
            num += 1
        
        #信息素挥发
        for i in range(0,len(self.customers)-1):
            for j in range(i+1,len(self.customers)):
                self.pheromone[i][j] *= (1 - self.sita)
                self.pheromone[j][i] = self.pheromone[i][j]

        if self.now_best.totalcost < self.best_solution.totalcost:
            total_distance = 0
            for route in self.now_best.routes:
                total_distance += route.distance

            delta = self.Q / total_distance
            self.best_solution = copy.deepcopy(self.now_best)
        else:
            delta = 0

        best_solution_node = []
        best_solution_path = []
        for i in self.best_solution.routes:
            for k in i.customers[:-1]:
                best_solution_node.append(k)
        best_solution_node.append(0)
        
        for i in range(0,len(best_solution_node)-1):
            best_solution_path.append([best_solution_node[i],best_solution_node[i+1]])

        for i in ant_solution_path:
            for j in i:
                if j in best_solution_path:
                    # 有改动
                    self.pheromone[j[0]][j[1]] += delta #* self.graph[j[0]][j[1]]
                    self.pheromone[j[1]][j[0]] = self.pheromone[j[0]][j[1]]

        for i in range(0,len(self.customers)-1):
            for j in range(i+1,len(self.customers)):
                if self.pheromone[i][j] > delta_max:
                    delta_max = self.pheromone[i][j]
                elif self.pheromone[i][j] < delta_min:
                    self.pheromone[i][j] = delta_min
                self.pheromone[j][i] = self.pheromone[i][j]
    
class Route:

    def __init__(self):

        super().__init__()
        self.load = 0
        self.distance = 0
        self.time = 0
        self.customers = []

class Solution:

    def __init__(self):
    
        super().__init__()
        self.totalcost = 0
        self.routes = []
        self.c1 = 0
        self.c2 = 0
        self.c3 = 0
        self.c4 = 0
        self.c5 = 0

class Customer:

    def __init__(self):

        super().__init__()
        self.id = 0         #某customer的id
        self.location_x = 0 #某customer横坐标
        self.location_y = 0 #某customer纵坐标
        self.q_i = 0        #某customer需求货量
        self.e_i = 0        #某customer时间窗最早
        self.l_i = 0        #某customer时间窗最晚
        self.t_load = 0     #装卸时间
        self.t_i = 0          #到达时间

if __name__ == "__main__":

    # 有改动
    z = [
        [0, [71.60, 96.72], 0, [0.0, 17.0], 0],
        [1, [70.70, 98.86], 1500, [6.0, 8.0], 1 / 3],
        [2, [70.47, 100.73], 500, [7.5, 9], 1 / 6],
        [3, [69.09, 99.42], 1500, [6.0, 8.0], 1 / 3],
        [4, [68.75, 98.41], 1500, [6.5, 8.33], 1 / 3],
        [5, [71.67, 101.61], 2000, [6.67, 8.5], 5 / 12],
        [6, [69.14, 101.44], 2000, [7.0, 9.0], 5 / 12],
        [7, [67.98, 100.32], 1800, [7.33, 9.0], 11 / 30],
        [8, [70.21, 102.49], 1000, [7.5, 9.0], 1 / 4],
        [9, [67.91, 98.22], 1000, [7.0, 8.5], 1 / 4],
        [10, [66.67, 100.79], 1000, [7.0, 9.0], 1 / 4],
        [11, [67.42, 102.81], 1000, [7.5, 9.5], 1 / 4],
        [12, [69.22, 103.54], 500, [7.5, 9.0], 1 / 6],
        [13, [65.98, 102.38], 500, [7.5, 9.5], 1 / 6],
        [14, [73.00, 101.03], 1500, [7.5, 9.0], 1 / 3],
        [15, [72.98, 102.44], 2000, [6.83, 8.5], 5 / 12],
        [16, [71.86, 103.30], 1500, [7.0, 8.67], 1 / 3],
        [17, [71.00, 102.40], 1500, [7.0, 8.67], 1 / 3],
        [18, [72.03, 101.11], 500, [7.83, 9.0], 1 / 6],
        [19, [69.82, 98.65], 2500, [6.5, 8.5], 1 / 2],
        [20, [71.21, 98.11], 1000, [7.83, 9.0], 1 / 4]
    ]
    hybrid_ant_colony = HybridAntColony(z)
    hybrid_ant_colony.init()
    hybrid_ant_colony.reset()
    hybrid_ant_colony.a_star()

    hybrid_ant_colony.init_phe()

    temp_solution = []
    convergence_times = 0
    convergence = {}
    convergence[convergence_times] = hybrid_ant_colony.best_solution.totalcost
    delta_max = hybrid_ant_colony.lambda_phe * hybrid_ant_colony.Q / (len(hybrid_ant_colony.customers) * (len(hybrid_ant_colony.customers)-1))
    delta_min = delta_max / 5
    for i in range(0,hybrid_ant_colony.ant):
        hybrid_ant_colony.reset()
        hybrid_ant_colony.construct_solution()
        temp_solution.append(hybrid_ant_colony.solution)
    last_best_totalcost = hybrid_ant_colony.best_solution.totalcost
    hybrid_ant_colony.update_pheromone(temp_solution,delta_max,delta_min)
    if hybrid_ant_colony.now_best.totalcost < last_best_totalcost:
        convergence_times = 1
        convergence[convergence_times] = hybrid_ant_colony.now_best.totalcost
        # distance = 0
        # for route in hybrid_ant_colony.now_best.routes:
        #     distance += route.distance
        # delta_max = (len(hybrid_ant_colony.customers)/20) / (1 - hybrid_ant_colony.sita) * 1 / distance
        # delta_min = delta_max / 5
    last_best_solution = copy.deepcopy(hybrid_ant_colony.now_best)

    num = 0
    for i in range(0,hybrid_ant_colony.times-1):

        if (i + 2) > hybrid_ant_colony.times / 5 and (i + 2) <= hybrid_ant_colony.times / 3:
            hybrid_ant_colony.sita = 0.5
        elif (i + 2) > hybrid_ant_colony.times / 3 and (i + 2) <= hybrid_ant_colony.times / 2:
            hybrid_ant_colony.sita = 0.3
        elif (i + 2) > hybrid_ant_colony.times / 2:
            hybrid_ant_colony.sita = 0.1

        temp_solution = []
        for j in range(0,hybrid_ant_colony.ant):
            hybrid_ant_colony.reset()
            hybrid_ant_colony.construct_solution()
            temp_solution.append(hybrid_ant_colony.solution)
        last_best_totalcost = hybrid_ant_colony.best_solution.totalcost
        hybrid_ant_colony.update_pheromone(temp_solution,delta_max,delta_min)
        if hybrid_ant_colony.now_best.totalcost < last_best_totalcost:
            convergence_times = i + 2
            convergence[convergence_times] = hybrid_ant_colony.now_best.totalcost
            # distance = 0
            # for route in hybrid_ant_colony.now_best.routes:
            #     distance += route.distance
            # delta_max = (len(hybrid_ant_colony.customers)/20) / (1 - hybrid_ant_colony.sita) * 1 / distance
            # delta_min = delta_max / 5
        
        if hybrid_ant_colony.now_best.totalcost != last_best_solution.totalcost:
            hybrid_ant_colony.q_0 *= (1 - (hybrid_ant_colony.now_best.totalcost - last_best_solution.totalcost) / last_best_solution.totalcost)
        else:
            if num < hybrid_ant_colony.N_MAX:
                num += 1
            else:
                hybrid_ant_colony.q_0 *= hybrid_ant_colony.eta
        
        last_best_solution = copy.deepcopy(hybrid_ant_colony.now_best)

    # 出路线图
    fig, ax = plt.subplots()
    Path = mpath.Path
    path_data = []

    for i in hybrid_ant_colony.best_solution.routes:
        for j in i.customers:
            if not path_data:
                path_data.append((Path.MOVETO,(hybrid_ant_colony.customers[j].location_x,hybrid_ant_colony.customers[j].location_y)))
            path_data.append((Path.LINETO,(hybrid_ant_colony.customers[j].location_x,hybrid_ant_colony.customers[j].location_y)))

    codes, verts = zip(*path_data)
    path = mpath.Path(verts, codes)
    x, y = zip(*path.vertices)
    line, = ax.plot(x, y, 'go-')
    ax.grid()
    ax.axis('equal')
    # plt.show()

    x_label = []
    y_label = []
    y_last_label = convergence[1]
    for i in convergence.keys():
        if i == 0 or i == 1:
            x_label.append(i)
            y_label.append(convergence[i])
        else:
            if (i-1) in convergence.keys():
                x_label.append(i)
                y_label.append(convergence[i])
                y_last_label = convergence[i]
            else:
                x_label.append(i-1)
                y_label.append(y_last_label)
                x_label.append(i)
                y_label.append(convergence[i])
                y_last_label = convergence[i]
    x_label.append(200)
    y_label.append(hybrid_ant_colony.best_solution.totalcost)
    plt.plot(x_label, y_label)
    plt.xlabel('Iterations')
    plt.ylabel('Cost')
    # plt.show()

    with open('C:\\Users\\98600\\Desktop\\new_hybrid_ant_colony.txt', 'a', encoding='utf-8') as fp:
        for i in hybrid_ant_colony.best_solution.routes:
            fp.write(str(i.customers)+' ')
        fp.write(str(hybrid_ant_colony.best_solution.totalcost)+' ')
        fp.write(str(hybrid_ant_colony.best_solution.c1)+' ')
        fp.write(str(hybrid_ant_colony.best_solution.c2)+' ')
        fp.write(str(hybrid_ant_colony.best_solution.c3)+' ')
        fp.write(str(hybrid_ant_colony.best_solution.c4)+' ')
        fp.write(str(hybrid_ant_colony.best_solution.c5)+' ')
        fp.write(str(convergence_times)+'\n')