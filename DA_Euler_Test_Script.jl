using NLsolve
using Plots
gr()

include("DAEs.jl")
include("Discrete_Methods.jl")

#START WITH ROBERTS
#Define Values for roberts
#Define Params and initial values
x0 = [1.,0.,0.]
params = 0. #Placeholder
t0exp = -6. #Base ten exponent of starting time
tfexp = 5. #Base ten exponent of ending time
steps = 10000 #number of Euler timesteps
EulerStepsPerPlotStep = 100
timesteps = 10 .^(collect(range(t0exp,stop=tfexp,length=steps)))#prep a vector of timesteps



#solve roberts with Discrete Euler
#x =discreteEuler(DAEroberts!,params,x0,timesteps)

#plot the solution
#plot(timesteps[1:25:end],x[:,1:25:end]', xscale =:log10,xlims =(.000001,100000),layout=(3,1))
#savefig("Robertson.png")

#RUN PENDULUM TEST
mass1 = 3.
mass2 = 1.
L1 = 1.
L2 = 3.
g = 9.81
tmax = 60
stepspersecond = 1000 #number of Euler timesteps
plotpointspersecond = 20
params = [mass1,mass2,L1,L2,g]
x0=[0,0,-9.,5.]

steps = stepspersecond*tmax
timesteps = collect(range(0,stop=tmax,length=steps))

x =discreteEuler(simplependulum!,params,x0,timesteps)

println("Simulation Complete")

#plot(x[3,1:EulerStepsPerPlotStep:steps],x[4,1:EulerStepsPerPlotStep:steps], xlims =(-2,2), ylims =(-2,2))
#plot!(@. sin(x[1,1:EulerStepsPerPlotStep:steps]),-cos(x[1,1:EulerStepsPerPlotStep:steps]))

#Create an animation
EulerStepsPerPlotStep = Int(ceil(stepspersecond/plotpointspersecond))
swing = Animation()
framerate =20
frames = Int(ceil(framerate * tmax))
timeperframe = 1/framerate
stepsperframe = Int(floor(stepspersecond*timeperframe))

for n = 1:1:frames
    th1 =x[1,1:EulerStepsPerPlotStep:n*stepsperframe]
    th2 =x[2,1:EulerStepsPerPlotStep:n*stepsperframe]
    x1 = L1*sin.(th1)
    y1 = -L1*cos.(th1)
    x2 = L1*sin.(th1).+L2*sin.(th2)
    y2 = -L1*cos.(th1) .- L2*cos.(th2)
    plot( last(x1,30),last(y1,30), xlims =(-(L1+L2),(L1+L2)), ylims =(-(L1+L2),(L1+L2)),size = (600,600))
    plot!([0,last(x1)],[0,last(y1)])
    plot!( last(x2,120),last(y2,120))
    plot!([last(x2),last(x1)],[last(y2),last(y1)])

    frame(swing)
end

println("Plot Complete")

gif(swing,"swing.gif",fps = framerate)
