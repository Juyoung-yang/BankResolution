# parameter values

qd = 0.9827
β = 0.9827
Rf = 0.0176 # net interest rate, whereas Rl is gross interest rate on loan 
wr = 0.428
α = 0.08
ρ = 0.4 # 은행 실패시의 패널티, calibration table의 η, 실패은행의 주가 하락률이 60%라고 가정했을 때, 1-η = 0.4
# 0.72
g = 0.1816246 # 예금부채 대비 보증자산 규모
ξ = 0.0214
cF = 0.0 # DELETE
dBar = 1.0
σ = 0.9932
τC = 0.215 # 최근 10년 실효 법인세 평균: 21.5%
z = 1.0

γ = 1.05 ## DELETE
ϕ = 1.0 ## DELETE
sigHat = 1.5; 
lconstr = 1.2;

δL = 27.5 # 24년도 예금규모 상위 25%: 27.5 조원 ## 원래는 10 
δM = 55.0 # 24년도 예금규모 상위 50%: 55.0 조원 ## 원래는 60
δH = 322.0 # 24년도 예금규모 상위 75%: 322.0 조원 ## 원래는 100

H = [0.947368421052632 0.0526315789473684 0; 0.0235294117647059 0.929411764705882 0.0470588235294118; 0 0.0945945945945946 0.905405405405405]  # transition matrix: [iDelta, iDeltaPrime]

# 고정이하여신비율 데이터를 통해 estimate

λL = 0.0047; # 2024년 고정이하여신비율 상위 33분위 
λM = 0.0074 # 2024년 고정이하여신비율 상위 66분위 
λH = 0.1796 # 평균 고정이하여신비율 if year == 1999 & bank that failed in early 2000 (서울,평화,조흥)

# F = [0.9997 0.0002 0.0001; 0.005 0.99 0.005; 0 0.1 0.9] # transition matrix: [iLambda, iLambdaPrime]
# F = [0.9 0.07 0.03; 0.07 0.9 0.03; 0 0.88 0.12] # transition matrix: [iLambda, iLambdaPrime]
F = [0.625 0.2875 0.0875; 0.22449 0.500 0.27551; 0.075949 0.341772 0.582279] # from 고정이하여신비율 

M = 1.0 # 정책변수: number of banks 

## calibrated parameter ##
cM = 1.3e-5
cO = 0.1 # 0.0 # no operation cost; # 0.2
cL = -0.003 # 0.1 # let operation cost for deposit being very small; 
# cL originally 0.5; with 0.5, optimization does not converge
ϵ = -1.1 # loan demand elasticity, - 1.1 for April
E = 100.0 # loan demand shifter

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

# calibrated parameters: 8+1 parameters
# F[1,3]
# F[2,3]
# F[3,3]
# cM
# cO
# cL
# ϵ
# E

# targeted moments: 8+1 moments
# 1. average loan-to-deposit ratio 
# 2. average debt-to-deposit ratio
# 3. average devidend-to-asset ratio
# 4. average loan interest revenue-to-loan ratio
# 5. eq. loan interest rate

# (regarding high loan default shock) moments conditional on bank failture?
# IMF 당시 파산한 은행의 balance sheet 
# 6. E(loan|lambda = λH & bank failure))
# 7. 
# 8.
# 9.