import csv
import xlwt
import numpy as np
import math
import random as rd
from scipy.stats import norm
from timeit import default_timer as timer

def ln(x):
    if x <= 0:
        return -float('inf')
    return math.log(x, math.exp(1))

def cumuDesity(x):
    return norm.cdf(x,0,1)

def getd1(S, K, T, sigma, mu):

    return (ln(float(S) / K) + (mu + 0.5 * sigma**2) * T) / (math.sqrt(T) * sigma)

def getd2(S, K, T, sigma, mu):

    return (ln(float(S) / K) + (mu + 0.5 * sigma**2) * T) / (math.sqrt(T) * sigma) - math.sqrt(T) * sigma

def geoAsianOption(S, sigma, r, T, K, N,option):
    K = float(K)
    N = float(N)
    new_sigma = sigma * math.sqrt((N + 1) * (2 * N + 1) / (6 * N**2))
    new_mu = (r - 0.5 * sigma**2) * (N + 1) / (2 * N) + 0.5 * new_sigma**2

    d1 = getd1(S, K, T, new_sigma, new_mu)
    d2 = getd2(S, K, T, new_sigma, new_mu)

    if option == 'C' or option == 'Call':
        return math.exp(-r * T) * (S * math.exp(new_mu * T) * cumuDesity(d1) - K * cumuDesity(d2))
    elif option == 'P' or option == 'Put':
        return math.exp(-r * T) * (- S * math.exp(new_mu * T) * cumuDesity(- d1) + K * cumuDesity(- d2))
    else:
        print "Please input correct option type: 'Call' or 'Put' "


def geoBasketOption(S1, S2, sigma1, sigma2, r, T, K, cor, option):
    sigma = math.sqrt(sigma1**2 + 2 * cor * sigma1 * sigma2 + sigma2**2) / 2
    mu = r - 0.5 * (sigma1**2 + sigma2**2) / 2 + 0.5 * sigma**2
    Bg = math.sqrt(S1 * S2)

    d1 = getd1(Bg, K, T, sigma, mu)
    d2 = getd2(Bg, K, T, sigma, mu)

    if option == 'C' or option == 'Call':
        return math.exp(-r * T) * (Bg * math.exp(mu * T) * cumuDesity(d1) - K * cumuDesity(d2))
    elif option == 'P' or option == 'Put':
        return math.exp(-r * T) * (- Bg * math.exp(mu * T) * cumuDesity(- d1) + K * cumuDesity(- d2))
    else:
        print "Please input correct option type: 'Call' or 'Put' "


def binoTreeMethod(S, sigma, r, T, K, N, option):
    dt = float(T) / N
    u = math.exp(sigma * math.sqrt(dt))
    d = 1 / u
    p = (math.exp(r * dt) - d) / (u - d)
    df = math.exp(-r * dt)

    st = [0.0] * (N + 1)   #stock price
    value = [0.0] * (N + 1)  # value of option
    for i in range(N + 1):
        st[i] = S * (u**(N - i)) * (d**(i))

    if option == 'C' or option == 'Call':
        for i in range(N + 1):
            value[i] = max(st[i] - K, 0)

        for j in range(N):  # iterate from the end
            for k in range(N - j):
                st[k] = st[k] / u
                value[k] = max(st[k] - K, df * (p * value[k] +(1 - p) * value[k + 1]))
        return value[0]

    if option == 'P' or option == 'Put':
        for i in range(N + 1):
            value[i] = max(- st[i] + K, 0)

        for j in range(N):  # iterate from the end
            for k in range(N - j):
                st[k] = st[k] / u
                value[k] = max(- st[k] + K, df * (p * value[k] +(1 - p) * value[k + 1]))
        return value[0]

def N (x):
    return (1/math.sqrt(2*math.pi))*(math.pow(math.e,-x*x/2))

# Block_Scholes
def cdf (x):
    if(x < -4.9):
        return 0
    else:
        if(x > 4.9):
            return 1
    f = 0.0
    pc = -5.0
    step = 0.00001
    while pc<x:
        f += N(pc) * step
        pc = pc + step
    return f

def Block_Scholes_call (s,k,T,v,repo,rfr):
    C=0.0
    t=0.0
    d1 = (math.log(s / k) + (rfr-repo) * (T - t)) / (v * math.sqrt(T - t)) + (0.5) * v * math.sqrt(T - t)
    d2 = (math.log(s / k) + (rfr-repo) * (T - t)) / (v * math.sqrt(T - t)) - (0.5) * v * math.sqrt(T - t)
    C = (s*math.exp((0-repo)*(T-t))*cdf(d1)) - (k*math.pow(math.e,(-rfr)*(T-t)))*cdf(d2)
    return C

def Block_Scholes_put (s,k,T,v,repo,rfr):
    P=0.0
    t = 0.0
    d1 = (math.log(s / k) + (rfr-repo) * (T - t)) / (v * math.sqrt(T - t)) + (0.5) * v * math.sqrt(T - t)
    d2 = (math.log(s / k) + (rfr-repo) * (T - t)) / (v * math.sqrt(T - t)) - (0.5) * v * math.sqrt(T - t)
    P = k * math.pow(math.e,-rfr*(T-t)) * cdf(-d2) - s * math.exp(-repo*(T-t))* cdf(-d1)
    return P

def vol_call(s,k,q,T,r,v):
    # s: underlying price, k:strike price, q:borrowing cost and dividend, T:time to expire
    # r: risk-free interest rate, v:market price
    s = float(s)
    k = float(k)
    q = float(q)
    T = float(T)
    r = float(r)
    v = float(v)
    vol = 2*abs(np.sqrt(np.log(s/k)+(r-q)*(T))/T) #inital volatility
    upper = 1
    lower =0
    C=0
    while abs(C-v)>1e-7:
        d1 = (np.log(s/k)+(r-q+vol**2/2)*T) / (vol*np.sqrt(T))
        d2 = d1 - vol*np.sqrt(T)
        C = s * np.exp(-q*T)*norm.cdf(d1) - k*np.exp(-r*T)*norm.cdf(d2)
        if C-v > 0:
            upper = vol
            vol = (vol+lower)/2
        else:
            lower = vol
            vol = (vol + upper)/2
    return vol

def vol_put(s,k,q,T,r,v):
    # s: underlying price, k:strike price, q:borrowing cost and dividend, T:time to expire
    # r: risk-free interest rate, v:market price
    s = float(s)
    k = float(k)
    q = float(q)
    T = float(T)
    r = float(r)
    v = float(v)
    vol = 2*abs(np.sqrt(np.log(s/k)+(r-q)*(T))/T) #inital volatility
    upper = 1
    lower =0
    P=0
    while abs(P-v)>1e-7:
        d1 = (np.log(s/k)+(r-q+vol**2/2)*T) / (vol*np.sqrt(T))
        d2 = d1 - vol*np.sqrt(T)
        P  =  k*np.exp(-r*T)*norm.cdf(-d2) - s * np.exp(-q*T)*norm.cdf(-d1)
        if P-v > 0:
            upper = vol
            vol = (vol+lower)/2
        else:
            lower = vol
            vol = (vol + upper)/2
    return vol

def BWM(S, dt, c1, c2, random):
    return S * np.exp(c1 * dt + c2 * random)


def arith_asian_option(S, sigma, r, T, K, step, type, path, cv):
    np.random.seed(0)
    paths = np.zeros((path, step), order='F')
    random = np.zeros((path, step), order='F')
    for i in range(0, path):
        random[i, :] = np.random.standard_normal(step)
    c1 = r - 0.5 * sigma * sigma
    dt = float(T) / step
    c2 = sigma * np.sqrt(dt)
    paths[:, 0] = BWM(S, dt, c1, c2, random[:, 0])
    for i in range(1, step):
        s = paths[:, i - 1]
        paths[:, i] = BWM(s, dt, c1, c2, random[:, i])

    arithMean = paths.mean(1)
    logPaths = np.log(paths)
    geoMean = np.exp(1 / float(step) * logPaths.sum(1))

    if type == 'C':
        arith_payoff = np.maximum(arithMean - K, 0)*np.exp(-r*T)
        geo_payoff = np.maximum(geoMean - K, 0)*np.exp(-r*T)
    elif type == 'P':
        arith_payoff = np.maximum(K - arithMean, 0)*np.exp(-r*T)
        geo_payoff = np.maximum(K - geoMean, 0)*np.exp(-r*T)
    else:
        return 404

    if cv == 'NULL':
        return np.mean(arith_payoff)

    else:
        XY = arith_payoff*geo_payoff
        covXY = np.mean(XY) - (np.mean(geo_payoff) * np.mean(arith_payoff))
        theta = covXY/np.var(geo_payoff)
        geo = geoAsianOption(S, sigma,r, T, K, step, type)
        z = arith_payoff + theta*(geo - geo_payoff)
        return np.mean(z)


def arith_basket(S1, S2, sigma1, sigma2, r, T, K, corr, type, path, cv):
    np.random.seed(0)
    z1 = np.random.standard_normal(path)
    z = np.random.standard_normal(path)
    z2 = corr*z1+math.sqrt(1-corr**2)*z
    S1_T = S1*np.exp((r-0.5*sigma1**2)*T+sigma1*np.sqrt(T)*z1)
    S2_T = S2*np.exp((r-0.5*sigma2**2)*T+sigma2*np.sqrt(T)*z2)
    ba_T = (S1_T+S2_T)/2
    bg_T = np.exp((np.log(S1_T)+np.log(S2_T))/2)
    if type == 'C':
        arith_payoff = (ba_T-K)*math.exp(-r*T)
        geo_payoff = (bg_T-K)*math.exp(-r*T)
    else:
        arith_payoff = (K-ba_T)*math.exp(-r*T)
        geo_payoff = (K-bg_T)*math.exp(-r*T)
    for i in range(0, path):
        arith_payoff[i] = max(arith_payoff[i], 0)
        geo_payoff[i] = max(geo_payoff[i], 0)


    if cv == 'NULL':
        p_mean = np.mean(arith_payoff)
        p_std = np.std(arith_payoff)
        return p_mean

    else:
        XY = [0.0]*path
        for i in range(0, path):
            XY[i] = arith_payoff[i]*geo_payoff[i]
        covXY = np.mean(XY) - np.mean(arith_payoff)*np.mean(geo_payoff)
        theta = covXY/np.var(geo_payoff)

        geo = geoBasketOption(S1, S2, sigma1, sigma2, r, T, K ,corr, type)
        Z = arith_payoff + theta * (geo - geo_payoff)
        z_mean = np.mean(Z)
        z_std = np.std(Z)
        return z_mean
