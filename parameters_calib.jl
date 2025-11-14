pwd()

# for simulation
bigT = 100
bigN = 10
bigJ = 100
trim = 10

# for boundary of root finding
# a = 1.0+eps()
b = 2.0-eps()

# for root finding with initial guess 
a = 1.04+eps()

sigHat = 1.15
lconstr = 2.4
E = 210.0


# under 최종 parameter: sigHat = 1.15/ lconstr = 2.4/ E = 210.0 with cL = -0.005
# target moment: need to have 5 moments; 3 for balance sheet rated, loan interest rate, relative size of dividend 
debt_to_liability = 0.122136 ## 2000년 ~ 2024년 은행간 연 평균 -------- model: 0.45036739269891835
loan_to_asset = 0.713694 ## 2000년 ~ 2024년 은행간 연 평균 -------- model: 0.6268244707025968
capital_to_deposit = 0.0913644 ## 2000년 ~ 2024년 은행간 연 평균 --------- model: 0.2833193942901851
loan_rate = 1.0473  ## 2024년 기준 ------- model: 1.0402875686320894
dividend_to_deposit = .0035879 ## 2024년 기준 : 0.35% ------ model: -0.031056317896547517 (전체)/ 0.0553845549300786 (배당지급 하는 경우에만)
## bank failure rate = 3.8% or 0.038 ------- model: 0.08379545454545455
## 예대율 비율: ----- model : 1.3753338636637973

