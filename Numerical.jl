using Distributed
#addprocs(4)
@everywhere using NLsolve
@everywhere using Parameters
@everywhere using LinearAlgebra
@everywhere using DifferentialEquations
@everywhere using ParameterizedFunctions
@everywhere using Plots
@everywhere using PyPlot
@everywhere using FileIO
@everywhere using CSV
@everywhere using DelimitedFiles
@everywhere using DataFrames
using Plots.PlotMeasures
using LaTeXStrings
#@everywhere #using MatPlotLib
@everywhere include("parameters.jl")


#### RUnge_kutta
function eqODE(du,u,P,t)
    du[1] = P.n*P.phi*(1-u[1]/P.K)*u[3]-P.gammaT*u[1]-P.rho*u[5]*u[1]-P.alphaE*u[6]*u[1]-P.muET*u[1] ##ET
    du[2] = P.gammaT*u[1]-P.alphaN*u[6]*u[2]-P.eta*u[2]-P.muNT*u[2] ##NT
    du[3] = P.eta*u[2]-P.muAT*u[3] ##AT
    du[4] = P.rho*u[5]*u[1]-P.gammaP*u[4]-P.alphaP*u[6]*u[4]-P.muEP*u[4] ##EP
    du[5] = P.gammaP*u[4]-P.muAP*u[5] ##AP
    du[6] = P.betaE*u[6]*u[1]+P.betaN*u[6]*u[2]+P.betaP*u[6]*u[4]-P.muAM*u[6] ##AM
end

#param = CSV.File("Parametros5000.csv",header = false,type = Float64)
param = readdlm("parametros600.dat",header = false)
sim = 10
P = EDOparameters(
    n = param[sim,2],
    phi = param[sim,3],
    K = 1,#param[sim,4],
    gammaT = param[sim,5],
    rho = param[sim,6],
    alphaE = param[sim,7],
    muET = param[sim,8],
    alphaN = param[sim,9],
    eta = param[sim,10],
    muNT = param[sim,11],
    muAT = param[sim,12],
    gammaP = param[sim,13],
    alphaP = param[sim,14],
    muEP = param[sim,15],
    muAP = param[sim,16],
    betaE = param[sim,17],
    betaN = param[sim,18],
    betaP = param[sim,19],
    muAM = param[sim,20]
)


u0=[2.5;0.00;0.00;0.0;0.0;3.0]
P = EDOparameters()

Q0 = (P.eta*P.phi*P.gammaT*P.n)/(P.muAT*(P.gammaT+P.muET)*(P.eta+P.muNT))
ET0 = P.K/Q0*(Q0-1)
NT0 = P.gammaT/(P.eta+P.muNT)*ET0
Qp = (P.rho*P.gammaP)/(P.muAP*(P.gammaP+P.muEP))*ET0
Qpre = (P.betaE/P.muAM*ET0+P.betaN/P.muAM*NT0)

#########


tspan=(0.0,1000.0)
prob=DifferentialEquations.ODEProblem(eqODE,u0,tspan,P)
sol=DifferentialEquations.solve(prob,saveat=0.1)
#sol[6,:]
sol2 = zeros(Float64,size(sol)[2],size(sol)[1])
for i = 1:size(sol)[1]
    sol2[:,i] = sol[i,:]
end

lb = ["Tuta eggs" "Tuta larvae" "Tuta adults" "Parasitoid eggs" "Parasitoid adults" "Predator"]
cl = [:black :grey :lightgrey :blue :red :green]

Plots.plot(sol.t,sol2[:,[6]],ylabel = "Number of individuals",xlabel = "Time (days)",title = "Q_0<1",labels = lb,lw = 3,c = cl,xlim = (0,1000))
sol2[end,:]


P.K = 0.003
P.n = 50.0
P.alphaE = 0.1
P.muAT = 1/16
P.muET = 0
P.muNT = 0
tspan=(0.0,1000.0)
prob=DifferentialEquations.ODEProblem(eqODE,u0,tspan,P)
sol=DifferentialEquations.solve(prob,saveat=0.1)
#sol[6,:]
sol2 = zeros(Float64,size(sol)[2],size(sol)[1])
for i = 1:size(sol)[1]
    sol2[:,i] = sol[i,:]
end

lb = ["Tuta eggs" "Tuta larvae" "Tuta adults" "Parasitoid eggs" "Parasitoid adults" "Predator"]
cl = [:black :grey :lightgrey :blue :red :green]

Plots.plot(sol.t,sol2[:,[4;6]],ylabel = "Number of individuals",xlabel = "Time (days)",title = "Q_0<1",labels = lb,lw = 3,c = cl,xlim = (0,1000))
sol2[end,:]

Plots.savefig("q_greater_one.png")
Plots.plot(sol.t,sol[5,:],yaxis = :log)

cl = [:black :grey :lightgrey]
lb = ["Parasitoid eggs" "Parasitoid adults" "Predator"]
Plots.plot(sol.t,sol2[:,4:6],ylabel = "Number of individuals",xlabel = "Time (days)",labels = lb,lw = 3,c = cl,xlim = (0,10000))

Plots.plot(sol.t,sol2[:,6],ylabel = "Number of individuals",xlabel = "Time (days)",labels = "Predator",lw = 3,c = :black,xlim = (0,10000))


########################################
######### Stability Regions#########
###################################


function eqODE(du,u,P,t)
    du[1] = P.n*P.phi*(1-u[1]/P.K)*u[3]-P.gammaT*u[1]-P.rho*u[5]*u[1]-P.alphaE*u[6]*u[1]-P.muET*u[1]
    du[2] = P.gammaT*u[1]-P.alphaN*u[6]*u[2]-P.eta*u[2]-P.muNT*u[2]
    du[3] = P.eta*u[2]-P.muAT*u[3]
    du[4] = P.rho*u[5]*u[1]-P.gammaP*u[4]-P.alphaP*u[6]*u[4]-P.muEP*u[4]
    du[5] = P.gammaP*u[4]-P.muAP*u[5]
    du[6] = P.p_t*P.alphaE*u[6]*u[1]+P.p_t*P.alphaN*u[6]*u[2]+P.p_pa*P.alphaP*u[6]*u[4]-P.muAM*u[6]
end

function run_several_times(pp,rho_v,time_run::Float64)

    qpre_v  = zeros(Float64,length(pp)*length(rho_v))
    qp_v = zeros(Float64,length(pp)*length(rho_v))
    point = zeros(Int64,length(pp)*length(rho_v))
    Et = zeros(Float64,length(pp)*length(rho_v))
    Nt = zeros(Float64,length(pp)*length(rho_v))
    At = zeros(Float64,length(pp)*length(rho_v))

    k::Int64 = 1
    for i = 1:length(rho_v)
        for j = 1:length(pp)

            println("$i $j")

            P = EDOparameters(
                rho = rho_v[i],
                p_t = pp[j],
                p_pa = 0.00125
            )


            u0=[20.0;10.0;10.0;1.0;0.0;1.0]

            Q0 = (P.eta*P.phi*P.gammaT*P.n)/(P.muAT*(P.gammaT+P.muET)*(P.eta+P.muNT))
            ET0 = P.K/Q0*(Q0-1)
            NT0 = P.gammaT/(P.eta+P.muNT)*ET0
            qp_v[k] = (P.rho*P.gammaP)/(P.muAP*(P.gammaP+P.muEP))*ET0
            qpre_v[k] = P.p_t*(P.alphaE/P.muAM*ET0+P.alphaN/P.muAM*NT0)

            tspan=(0.0,time_run)
            prob=DifferentialEquations.ODEProblem(eqODE,u0,tspan,P)
            sol=DifferentialEquations.solve(prob,saveat=0.1)

            if sol[4,end] < 0.00001 
                if sol[6,end] < 0.00001
                    point[k] = 1 ##only Tuta
                else
                    point[k] = 3 ## predator
                end
            elseif sol[6,end] < 0.00001
                point[k] = 2 ##parasitoid
            else
                point[k] = 4 #coexist
            end

            Et[k] = sol[1,end]
            Nt[k] = sol[2,end]
            At[k] = sol[3,end]

            k+=1

        end
    end

    return qp_v,qpre_v,point,Et,Nt,At

end

time_run = 10000.0


pp = 0:(0.000000365/3):0.00001825 #aproveitamento energetico
rho_v = 0:(0.0000213/3):0.001065
qp_v,qpre_v,point,Et,Nt,At = run_several_times(pp,rho_v,time_run)

writedlm("data_stability.dat",[qp_v qpre_v point Et Nt At])

cl = [:black; :blue; :red; :green]
cl_2 = cl[point]

Plots.scatter(qp_v,qpre_v,marker=(0.99,Plots.stroke(0)),c = cl_2,labels = "",xlabel = L"Q_{p}",ylabel=L"Q_{pre}")

Plots.savefig("stability_regions.png")

###drawing Line
rho_v = 0.0003124
qpre_v_line = zeros(Float64,length(rho_v))
qp_v_line = zeros(Float64,length(rho_v))
#q0_v_line = zeros(Float64,length(rho_v))
P = EDOparameters(p_pa = 0.00125)

for i = 1:length(rho_v)
    Q0 = (P.eta*P.phi*P.gammaT*P.n)/(P.muAT*(P.gammaT+P.muET)*(P.eta+P.muNT))
    ET0 = P.K*(Q0-1)/Q0
    qp_v_line[i] = (rho_v[i]*P.gammaP)*ET0/(P.muAP*(P.gammaP+P.muEP))
    ET1 = P.muAP*(P.gammaP+P.muEP)/(rho_v[i]*P.gammaP)
    Ep1 = (P.n*P.gammaT*P.eta*P.phi*P.muAP^2*(P.muEP+P.gammaP))*(qp_v_line[i]-1)/(P.K*rho_v[i]^2*P.gammaP^2*P.muAT*(P.eta+P.muNT))#Q0*(P.gammaT+P.muET)*P.muAP*P.muAP*(P.gammaP+P.muEP)*(qp_v_line[i]-1)/(P.K*rho_v[i]^2*P.gammaP^2)
    #NT0 = P.gammaT/(P.eta+P.muNT)*ET0
    qpre_v_line[i] = ET0/ET1*(1-P.p_pa*P.alphaP*Ep1/P.muAM)
end

pos = findall(x->x>=1.0,qp_v_line)
qp_v_line = qp_v_line[pos]
qpre_v_line = qpre_v_line[pos]

for i = 1:length(qpre_v_line)
    qpre_v_line[i] = qpre_v_line[i] > 0 ? qpre_v_line[i] : 0
end

pos = findall(x->x>=0.0,qpre_v_line)
qp_v_line = qp_v_line[pos]
qpre_v_line = qpre_v_line[pos]

Plots.plot!(qp_v_line,qpre_v_line,lw = 3,c =:black,labels = "")

Plots.plot!(qp_v_line,qp_v_line,lw = 3,c =:black,labels = "")

Plots.plot!(qp_v_line,ones(Float64,length(qp_v_line)),lw = 3,c =:black,labels = "")

Qpre = P.p*(P.alphaE/P.muAM*ET0+P.alphaN/P.muAM*NT0)



######################################
###Number of tuta#
#need to run the above code until "drawing line" starts

using Colors
using ColorSchemes
using ColorSchemeTools
using GR
using ColorBrewer
using Gnuplot
using Plots
using PyPlot
using DelimitedFiles
using PyCall
using FileIO

data = readdlm("data_stability.dat",header=false)

qp_v = data[:,1]
qpre_v = data[:,2]
point = data[:,3]
Et = data[:,4]
Nt = data[:,5]
At = data[:,6]

p1=Plots.scatter(qp_v,qpre_v,zcolor = At+Nt,m = (ColorSchemes.vik[:], 0.8, Plots.stroke(1, :black)),xlabel = L"Q_{p}",ylabel=L"Q_{pre}",label = "")
FileIO.save("number_of_tuta_t.png",p1)

i = 0:1:150
j = map(x->x*151+1,i)
x = qp_v[j]
y = qpre_v[1:151]

p_1 = 0:0.001:1
 
m_result = reshape(At+Nt,length(x),length(y))

gr()
p1=Plots.heatmap(x, y, m_result, aspect_ratio = 1,xlim=(0,5),c=cgrad(get(colorschemes[:hot],p_1)),xlab = L"Q_{p}",ylab = L"Q_{pre}")


Plots.plot!(qp_v_line,qpre_v_line,lw = 1.5,c =:black,labels = "")

Plots.plot!(qp_v_line,qp_v_line,lw = 1.5,c =:black,labels = "")

FileIO.save("number_of_tuta_heatmap_p000125.png",p1)


######## PHASE PORTRAIT ######
pyplot()

initial_conditions = zeros(Float64,7)
initial_conditions = [1.0;0.0;0.0;0.0;1.0;0.0;0.0]

for Yh = 0.0:0.1:0.9
    for Ys = 0.0:0.1:0.9
        global initial_conditions = [initial_conditions [1.0-Yh;Yh;0.0;0.0;1.0-Ys;Ys;0.0]]
    end
end
plotting_var = [2,6]
tspan=(0.0,300.0)

function run_solver(eqODE,u0,tspan,P)
    prob=DifferentialEquations.ODEProblem(eqODE,u0,tspan,P)
    sol=DifferentialEquations.solve(prob,saveat = 0.01)
    return sol
end

P = EDOparameters(beta1=0.80,gamma=500,step = 10000,alpha2=20)

sol = map(x->run_solver(eqODE,initial_conditions[:,x],tspan,P),1:size(initial_conditions)[2])

p=Plots.plot(sol[1][plotting_var[1],:],sol[1][plotting_var[2],:],c=:grey,label="",size = [1500,800])

for i = 2:size(initial_conditions)[2]
    global p=Plots.plot!(sol[i][plotting_var[1],:],sol[i][plotting_var[2],:],c=:grey,label="")
end
#p=Plots.xaxis!(scale=:log10)
#p=Plots.yaxis!(scale=:log)
p

function run_gradient(solucao,plotting_var,P,step_aux::Int64)
    dU_time_x = zeros(Float64,Int(floor(size(solucao)[1]/step_aux))+1)
    dU_time_y = zeros(Float64,Int(floor(size(solucao)[1]/step_aux))+1)
    U_time_x = zeros(Float64,Int(floor(size(solucao)[1]/step_aux))+1)
    U_time_y = zeros(Float64,Int(floor(size(solucao)[1]/step_aux))+1)
    dU = zeros(Float64,7)
    u = zeros(Float64,7)
    count::Int64 = 1

    for i = 1:step_aux:Int(size(solucao)[1])
        println(Int(i))
        u = solucao[Int(i)][:]
        U_time_x[count] = u[plotting_var[1]]
        U_time_y[count] = u[plotting_var[2]]
        eqODE(dU,u,P,0)
        dU_time_x[count] = dU[plotting_var[1]]/norm(dU)
        dU_time_y[count] = dU[plotting_var[2]]/norm(dU)
        count+=1
    end

    return U_time_x,U_time_y,dU_time_x,dU_time_y
end

size_quiver = 0.0005

for i = 1:size(initial_conditions)[2]

    U_time_x,U_time_y,dU_time_x,dU_time_y = run_gradient(sol[i].u,plotting_var,P,P.step)
    z = dU_time_x.*dU_time_x+dU_time_y.*dU_time_y
    p=Plots.quiver!(U_time_x,U_time_y,quiver=(dU_time_x*size_quiver,dU_time_y*size_quiver),label="",c = :gray)
end
p

p = Plots.yaxis!(ylim=(0,0.4))

p = Plots.xaxis!(xlim=(0,0.5))


file_index = "beta-1-case-3-gamma-menor-t"

FileIO.save(string("Phase_Portrait-","$(file_index)",".pdf"),p)


###bacias de atracao
@everywhere function eqODE(du,u,P,t)
   
    du[1] = P.Ah-P.beta1*u[6]*u[1]-P.muH*u[1]+P.muW1*u[2]
    du[2] = P.beta1*u[6]*u[1]-P.beta2*u[6]*u[2]-P.muH*u[2]-P.muW1*u[2]+P.muW2*u[3]
    du[3] = P.beta2*u[6]*u[2]-P.beta3*u[6]*u[3]-P.muH*u[3]-P.muW2*u[3]+P.muW3*u[4]
    du[4] = P.beta3*u[6]*u[3]-P.muH*u[4]-P.muW3*u[4]
    du[5] = P.As-P.betaS*u[7]*u[5]-P.muS*u[5]
    du[6] = P.betaS*u[7]*u[5]-P.muS*u[6]
    du[7] = P.alpha1*u[2]+P.alpha2*u[3]+P.alpha3*u[4]+P.gamma*u[6]*u[4]-P.muM*u[7]
end


@everywhere function run_points(x,initial_conditions,threshold_ys,P::EDOparameters)
    println("$x")

    tspan=(0.0,1000.0)
    prob=DifferentialEquations.ODEProblem(eqODE,initial_conditions,tspan,P)
    sol=DifferentialEquations.solve(prob,saveat = 0.1)

    e_point=sol.u[end][6]
    convergence_::Int64 = 0
    
    if e_point < threshold_ys
        convergence_ = 1
    else 
        convergence_ = 2
    end
    return convergence_,e_point
end

P = EDOparameters(beta1=0.7,gamma=650,alpha2=25)

initial_conditions = zeros(Float64,7)
initial_conditions = [0.203213;0.9;0.277217;0.0819154;0.639077;0.9;0.0225902]
initial_conditions2 = [0.203213;0.9;0.277217;0.0819154;0.639077;0.9;0.0225902]

plotting_var = [2,6]

for Yh = 0.0:0.005:0.9
    for Ys = 0.0:0.005:0.9
        initial_conditions2[plotting_var[1]] = Yh
        initial_conditions2[plotting_var[2]] = Ys
        global initial_conditions = [initial_conditions initial_conditions2]
    end
end


threshold_ys = 0.360923
result_ = pmap(x->run_points(x,initial_conditions[:,x],threshold_ys,P),1:size(initial_conditions)[2])
convergences_ = map(x->result_[x][1],1:size(result_)[1])

color_ = [:lightgrey :grey]

Plots.scatter(initial_conditions[plotting_var[1],:],initial_conditions[plotting_var[2],:],c = color_[convergences_],label = "",marker=(0.99,Plots.stroke(0)),size = [1500,800],bottom_margin =20mm,left_margin=20mm)

initial_conditions2 = [0.203213;0.9;0.277217;0.0819154;0.639077;0.9;0.0225902]
tspan=(0.0,1000)
initial_=[0.6;0.8]
initial_conditions2[plotting_var[1]] = initial_[1]
initial_conditions2[plotting_var[2]] = initial_[2]
prob=DifferentialEquations.ODEProblem(eqODE,initial_conditions2,tspan,P)
sol=DifferentialEquations.solve(prob,saveat = 0.1)

Plots.plot!(sol[plotting_var[1],:],sol[plotting_var[2],:],width=3,c=:black,label="")

p1 = xaxis!(xlabel = L"Y_{h1}",xguidefontsize=30)
p1 = yaxis!(ylabel = L"Y_s",yguidefontsize=30,yguidefonthalign = "left")


FileIO.save(string("Bacias-atracao","case-3-alpha-2-25-file-2",".png"),p1)




#########################################################################
#=   Figure for revision =#
    #######################


@everywhere include("parameters.jl")


#### RUnge_kutta
function eqODE(du,u,P,t)
    du[1] = P.n*P.phi*(1-u[1]/P.K)*u[3]-P.gammaT*u[1]-P.rho*u[5]*u[1]-P.alphaE*u[6]*u[1]-P.muET*u[1] ##ET
    du[2] = P.gammaT*u[1]-P.alphaN*u[6]*u[2]-P.eta*u[2]-P.muNT*u[2] ##NT
    du[3] = P.eta*u[2]-P.muAT*u[3] ##AT
    du[4] = P.rho*u[5]*u[1]-P.gammaP*u[4]-P.alphaP*u[6]*u[4]-P.muEP*u[4] ##EP
    du[5] = P.gammaP*u[4]-P.muAP*u[5] ##AP
    du[6] = P.betaE*u[6]*u[1]+P.betaN*u[6]*u[2]+P.betaP*u[6]*u[4]-P.muAM*u[6] ##AM
end

P = EDOparameters()

u01=[2.5;0.00;0.00;10.0;0.0;0.0] ##parasitoide
u02=[2.5;0.00;0.00;0.0;0.0;3.0]  ### predador
u03=[2.5;0.00;0.00;10.0;0.0;3.0] ## dois

##QP = 1.05, P.rho = 0.000223583  / Qp = 3.7 => P.rho = 0.000787863
##Qpre = 1.05, P.p_t = 0.00000383 / Qpre = 2.8 = > P.p_t = 0.000010213

#param = [[0.000223583;0.000010213],[0.000787863;0.00000383],[0.000787863;0.000010213]]


##QP = 2.0, P.rho = 0.000425872  / Qp = 3.7 => P.rho = 0.000787863
##Qpre = 0.575, P.p_t = 0.000001915 / Qpre = 2.8 = > P.p_t = 0.000010213

param = [[0.000425872;0.000001915],[0.000787863;0.000001],[0.000787863;0.000010213]]

index = 3

P.p_t = param[index][2]
P.rho = param[index][1]
P.betaE = P.p_t*P.alphaE
P.betaN = P.p_t*P.alphaN
P.betaP = P.p_pa*P.alphaP


Q0 = (P.eta*P.phi*P.gammaT*P.n)/(P.muAT*(P.gammaT+P.muET)*(P.eta+P.muNT))
ET0 = P.K/Q0*(Q0-1)
NT0 = P.gammaT/(P.eta+P.muNT)*ET0
Qp = (P.rho*P.gammaP)/(P.muAP*(P.gammaP+P.muEP))*ET0
Qpre = (P.betaE/P.muAM*ET0+P.betaN/P.muAM*NT0)

ETP = P.muAP*(P.gammaP+P.muEP)/(P.rho*P.gammaP)
NTP = P.gammaT/(P.eta+P.muNT)*ETP
EPP = P.n*P.phi*P.eta*P.muAP/(P.K*P.rho*P.gammaP*P.muAT)*NTP*(Qp-1)

Qall = (P.betaE/P.muAM*ETP+P.betaN/P.muAM*NTP+P.betaP/P.muAM*EPP)

P.eta/P.muAT



tspan=(0.0,10000.0)
#########
prob=DifferentialEquations.ODEProblem(eqODE,u01,tspan,P)
sol=DifferentialEquations.solve(prob,saveat=0.1)
#sol[6,:]
a = Tables.table(hcat(sol.u...)')
df = DataFrame(a)
colnames = [:ET;:NT;:AT;:EP;:AP;:AM]
rename!(df,colnames)
df.time = sol.t
Plots.plot(df.time,df.AP)

CSV.write("response/solution_$(index)_3.csv",df)

Plots.plot(df.AM,df.ET)



#########################################################################
#=   Hopf =#
    #######################
    
    using NLsolve
    using Parameters
    using LinearAlgebra
    using DifferentialEquations
    using ParameterizedFunctions
    using Plots
    using PyPlot
    using FileIO
    using CSV
    using DelimitedFiles
    using DataFrames
    using Plots.PlotMeasures
    using LaTeXStrings

    include("parameters.jl")
    include("equations.jl")

    
    P = EDOparameters()

    u03=[2.5;0.00;0.00;10.0;0.0;3.0] ## dois
    
    ##QP = 1.05, P.rho = 0.000223583  / Qp = 3.7 => P.rho = 0.000787863
    ##Qpre = 1.05, P.p_t = 0.00000383 / Qpre = 2.8 = > P.p_t = 0.000010213
    
    param = [[0.000223583;0.000010213],[0.000787863;0.00000383],[0.000787863;0.000010213],[0.000787863;0.000007295]]
    
    index = 4
    

    aux = 0.1:0.04:0.26
    m_a = zeros(Float64,1000001,0)
    m_b = zeros(Float64,1000001,0)
    for ind in 1:length(aux)

        P.p_t = param[index][2]
        P.rho = param[index][1]
        P.p_pa = aux[ind]

        P.betaE = P.p_t*P.alphaE
        P.betaN = P.p_t*P.alphaN
        P.betaP = P.p_pa*P.alphaP
        Q0 = (P.eta*P.phi*P.gammaT*P.n)/(P.muAT*(P.gammaT+P.muET)*(P.eta+P.muNT))
        ET0 = P.K/Q0*(Q0-1)
        NT0 = P.gammaT/(P.eta+P.muNT)*ET0
        Qp = (P.rho*P.gammaP)/(P.muAP*(P.gammaP+P.muEP))*ET0
        Qpre = (P.betaE/P.muAM*ET0+P.betaN/P.muAM*NT0)
        tspan=(0.0,1000000.0)

        #########
        prob=DifferentialEquations.ODEProblem(eqODE,u03,tspan,P)
        sol=DifferentialEquations.solve(prob,saveat=1.0)
        #sol[6,:]
        df = DataFrame(sol')
        colnames = [:ET;:NT;:AT;:EP;:AP;:AM]
        rename!(df,colnames)
        df.time = sol.t
       
        m_a = [m_a df.ET]
        m_b = [m_b df.AM]
    end


    #Plots.plot(m_a,m_b)

    Plots.plot(m_a[:,1],xlims = [0,250000],ylims = [80,150])
    #Plots.yaxis!(lims = [50;150])

    writedlm("hopf/solution_ET.dat",m_a)


    writedlm("hopf/solution_AM.dat",m_b)
    
    [i for i in aux]