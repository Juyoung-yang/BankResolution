# parameter values

qd = 0.9748
β = 0.9748
q_delta = 0.9801
Rf = 0.0259 # net interest rate, whereas Rl is gross interest rate on loan  2.59
cL = -0.005 # let operation cost for deposit being very small; 
cF = 1.0 # 필요없지만 일단 빼진 말자...
wr = 0.428
α = 0.08
ρ = 0.4 # 은행 실패시의 패널티, calibration table의 η, 실패은행의 주가 하락률이 60%라고 가정했을 때, 1-η = 0.4, 본문에서는 eta.
ρ_bailout = 0.72 # needs when calculating government spending on bailout
g = 0.1816246 # 예금부채 대비 보증자산 규모
ξ = 0.0214
τC = 0.205 # 최근 10년 실효 법인세 평균: 21.5%
lconstr = 1.2;
 
δL = 27.5 # 24년도 예금규모 상위 25%: 27.5 조원 ## 원래는 10 
δM = 55.0 # 24년도 예금규모 상위 50%: 55.0 조원 ## 원래는 60
δH = 322.0 # 24년도 예금규모 상위 75%: 322.0 조원 ## 원래는 100

H = [0.947368421052632 0.0526315789473684 0; 0.0235294117647059 0.929411764705882 0.0470588235294118; 0 0.0945945945945946 0.905405405405405]  # transition matrix: [iDelta, iDeltaPrime]

# 고정이하여신비율 데이터를 통해 estimate

λL = 0.0047; # 2024년 고정이하여신비율 상위 33분위 
λM = 0.0074 # 2024년 고정이하여신비율 상위 66분위 
λH = 0.1796 # 평균 고정이하여신비율 if year == 1999 & bank that failed in early 2000 (서울,평화,조흥)

dBar = 1.0
σ = 0.9932
cM = 1.3e-5
cO = 0.1
ϵ = -1.1 # (KEEP IT FIXED) loan demand elasticity, - 1.1 for April 

# alternative lambda values using 2024 data
# λL = 0.00275; # 2024년 고정이하여신비율 상위 25%
# λM = 0.00530 # 2024년 고정이하여신비율 상위 50%
# λH = 0.00745 # 2024년 고정이하여신비율 상위 75%

# F = [0.9997 0.0002 0.0001; 0.005 0.99 0.005; 0 0.1 0.9] # transition matrix: [iLambda, iLambdaPrime]
# F = [0.9 0.07 0.03; 0.07 0.9 0.03; 0 0.88 0.12] # transition matrix: [iLambda, iLambdaPrime]
F = [0.625 0.2875 0.0875; 0.22449 0.500 0.27551; 0.075949 0.341772 0.582279] # from 고정이하여신비율 
# F = [0.625 0.2875 0.0875; 0.22449 0.500 0.27551; 0.0 0.88 0.12] # from 고정이하여신비율 + April (for 3rd row)


M = 1.0 # 정책변수: number of banks 

###### main parameters to calibrate ######
E = 133.5 # loan demand shifter
sigHat = 1.5;


n_start = -10.0 # retained earning could be negative, in particular when bank doesn't fail and tax is too high 
n_npts = 5 # 11
n_stop = 100.0
l_start = eps()
l_npts = 5 # 20
l_stop = 100.0 # 200.0으로 높게 설정하면 boundary condition violated 
s_start = eps()
s_npts = 3 # 20
s_stop = 100.0
b_start = eps()
b_npts = 5 # 20
b_stop = 100.0

maxiter = 1000 # maximum number of iterations for every iteration
tol = 1e-2 # tolerance level for every iteration

