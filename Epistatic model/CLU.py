import numpy as np
import pandas as pd
import torch
import torch.optim as optim
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
from scipy.stats import iqr
import json
import gc
import warnings
from tqdm import tqdm  # 导入tqdm库

warnings.filterwarnings("ignore")
# 检查是否有可用的GPU
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# 读取数据并转换为Pandas数据框
data1_df = pd.read_csv('D:/一元抗寒/clean_data/ch.csv', index_col=0)
# 将Pandas数据框转换为Numpy数组
data1_np = data1_df.to_numpy()

# 将Numpy数组转换为PyTorch张量并移动到GPU
data1 = torch.tensor(data1_np, dtype=torch.float32).to(device)


# 辅助函数：计算多元正态分布的对数密度
def mahalanobis(x, center, cov):
    x_cen = x - center
    inv_cov = torch.inverse(cov)
    return torch.sum(x_cen @ inv_cov * x_cen, axis=1)


def dmvnorm(x, mean, sigma, log=True):
    distval = mahalanobis(x, mean, sigma)
    logdet = torch.sum(torch.log(torch.real(torch.linalg.eigvals(sigma))))
    log2pi = torch.tensor(1.8378770664093454835606594728112352797227949472755668)
    logretval = -((x.shape[1] * log2pi + logdet + distval) / 2)

    if log:
        return logretval
    else:
        return torch.exp(logretval)


# 辅助函数：线性方程计算
def linear_equation(x, linear_par):
    result = linear_par[:, 0][:, None] + linear_par[:, 1][:, None] * x
    return result


def logsumexp(v):
    vm = torch.max(v)
    return torch.log(torch.sum(torch.exp(v - vm))) + vm


def linear_equation_base(x, y):
    x = np.array(x, dtype=float).reshape(-1, 1)
    y = np.array(y, dtype=float)

    model = LinearRegression().fit(x, y)
    intercept = model.intercept_
    slope = model.coef_[0]

    return intercept, slope


def linear_equation_all(x, y, maxit=100):
    result = linear_equation_base(x, y)
    return result


# 辅助函数：计算协方差矩阵
def get_SAD1_covmatrix(par, n):
    phi = par[0]
    gamma = par[1]

    # 创建索引矩阵
    indices = torch.arange(1, n + 1, dtype=torch.float32).to(device)

    # 计算对角线元素
    diag_elements = (1 - torch.pow(phi, 2 * indices)) / (1 - torch.pow(phi, 2))

    # 计算非对角线元素
    sigma = torch.zeros((n, n), dtype=torch.float32).to(device)
    for i in range(n):
        sigma[i, i:] = torch.pow(phi, torch.arange(n - i, dtype=torch.float32).to(device)) * diag_elements[i]

    # 使矩阵对称
    sigma = sigma + sigma.T - torch.diag(torch.diag(sigma))

    # 缩放矩阵
    sigma *= torch.pow(gamma, 2)

    return sigma




def get_par_int(X, k, times1, n1):
    n, d = X.shape
    X1 = X

    cov_int = [0.9, iqr(np.diag(np.cov(X1.cpu().numpy(), rowvar=False)))]

    init_cluster = KMeans(n_clusters=k).fit(X.cpu())
    prob = np.bincount(init_cluster.labels_) / n
    times1 = times1.cpu().numpy()
  
    fit1 = np.array([linear_equation_all(times1, init_cluster.cluster_centers_[c, :n1]) for c in range(k)]).T
   
    return_obj = {
        'initial_cov_params': cov_int,
        'initial_mu_params': fit1.flatten(),
        'initial_probibality': prob
    }

    return return_obj


# 主函数
def Q_function(par, prob_log, omega_log, X, k, n1, times1):
    n = X.shape[0]
    n1, k = map(int, [n1, k])

    # 分割X矩阵
    X1 = X

    par_mu = par[2:]

    cov1 = get_SAD1_covmatrix(par[:2], n1)

    mu1_mat = torch.column_stack((par_mu[:k], par_mu[k:2 * k]))

    mu1 = linear_equation(times1, mu1_mat)
 
    mvn_log1 = torch.stack([dmvnorm(X1, mu1[i], cov1, True) for i in range(k)], dim=1)

    mvn_log = mvn_log1
    tmp = prob_log.T + mvn_log - omega_log
    Q = -torch.sum(tmp * torch.exp(omega_log))

    return Q


def fun_clu(data, k, Time1=None, trans=torch.log, initial_pars=None,
               iter_max=100):
    data = data[:, torch.argsort(torch.sum(data, axis=0))]
    print(data)
    n, d = data.shape
    print(n, d)
    n1 = data.shape[1]
    print(n1)
    iter = 0    
    epsilon = 10000
    # Time1 = Time1 if Time1 is not None else torch.arange(1, n1 + 1, dtype=torch.float32).to(device)  
    times1 = trans(torch.sum(data, axis=0) + 1).to(device)
    data = trans(data + 1)
    X1 = data

    
    print(times1)

    if initial_pars is None:
        par_int = get_par_int(data, k, times1, n1)
        prob_int = par_int['initial_probibality']
        cov_int = par_int['initial_cov_params']
        mu_int = par_int['initial_mu_params']
        initial_pars = np.hstack((cov_int, mu_int))
    else:
        prob_int = torch.tensor(initial_pars[:k]).to(device)
        initial_pars = initial_pars[k:]

    par = torch.tensor(initial_pars, dtype=torch.float32).to(device)
    par.requires_grad = True
    prob_logi = torch.log(torch.tensor(prob_int, dtype=torch.float32).to(device))

    while abs(epsilon) > 1000 and iter < iter_max:
        par_mui = par[2:]
        cov1i = get_SAD1_covmatrix(par[:2], n1).to(device).detach()
    
        mu1_mati = torch.column_stack((par_mui[:k], par_mui[k:2 * k])).to(device).detach()

        mu1i = linear_equation(times1, mu1_mati).to(device).detach()

        mvn_log1i = torch.stack([dmvnorm(X1, mu1i[i], cov1i, True) for i in range(k)], dim=1).to(device).detach()
        
        mvn_logi = mvn_log1i
        mvni = mvn_logi + prob_logi
        omega_logi = torch.stack([mvni[i] - torch.logsumexp(mvni[i], dim=0) for i in range(n)])
        
        LL_mem = Q_function(par, prob_logi, omega_logi, data, k, n1, times1)
        print(LL_mem)
        prob_expi = torch.stack([logsumexp(omega_logi[:, i]) for i in range(omega_logi.size(1))])
        prob_logi = prob_expi - torch.log(torch.tensor(n, dtype=torch.float32))
        optimizer = optim.AdamW([par], fused=True)

        for _ in tqdm(range(25), desc="Training Progress"):  # 为训练循环添加进度条
            optimizer.zero_grad()
            Q = Q_function(par, prob_logi, omega_logi, data, k, n1, times1)
            Q.backward()
            optimizer.step()
   
        par_hat = par
        par = par_hat
        LL_next = Q_function(par, prob_logi, omega_logi, data, k, n1, times1)
        epsilon = LL_next - LL_mem
        LL_mem = LL_next
        iter += 1
            
        print(f"Iter: {iter}, Log-Likelihood: {LL_next}")

    prob_log = prob_logi.detach().cpu().numpy()
    par = par.detach().cpu().numpy()
    LL_next = LL_next.detach().cpu().numpy()
    BIC = 2 * (LL_next) + 2*np.log(n) * (len(par) + k - 1)
    return {
        'cluster_number': k,
        'log-likelihood': LL_next,
        'BIC': BIC,
        'par': par,
        'prob_log': prob_log
    }

# result = terfun_clu(data1, data2, data3, k=95)
# result['par'] = result['par'].tolist()
# result['log-likelihood'] = result['log-likelihood'].tolist()
# result['prob_log'] = result['prob_log'].tolist()
#     # 动态生成保存路径
# save_path = f'/root/autodl-tmp/k95.json'
# with open(save_path, 'w') as f:
#     json.dump(result, f)

for k in range(5, 100,5):
    result = fun_clu(data1, k=k)
    result['par'] = result['par'].tolist()
    result['log-likelihood'] = result['log-likelihood'].tolist()
    result['prob_log'] = result['prob_log'].tolist()
    # 动态生成保存路径
    save_path = f'D:/一元抗寒/一元聚类/h_k{k}.json'
    with open(save_path, 'w') as f:
        json.dump(result, f)


# k_values = [5, 10, 15, 20, 25, 30, 45, 47]  # Replace with your desired values

# for k in k_values:
#     result = terfun_clu(data1, data2, data3, k=k)
#     result['par'] = result['par'].tolist()
#     result['log-likelihood'] = result['log-likelihood'].tolist()
#     result['prob_log'] = result['prob_log'].tolist()
    
#     # Dynamically generate the save path
#     save_path = f'/root/autodl-tmp/k{k}.json'
#     with open(save_path, 'w') as f:
#         json.dump(result, f)


