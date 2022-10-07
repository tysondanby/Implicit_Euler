#-----------DISCRETE EULER
function discreteEuler(residualfunction,parameters,xo,tsteps)#pass in the residual function, any known params, initial state (x0), and the points in time to evaluate at (tsteps)
    #outputs x, a vector of states (i.e. a vector of vectors), with index-by-index correspondence to tsteps.

    x = zeros(length(xo),length(tsteps)) #Allocate the x vector of state vectors
    x[:,1]=xo #X starts at x0

    xdotguess = zeros(length(xo)) #First guess of state rates for NLsolve to be used for every itteration.
    for n = 1:1:(length(tsteps)-1)#keep track of index
        time = tsteps[n]#get the time at the current step
        dt = tsteps[n+1]-tsteps[n] #get the change in time

        #Simplify the residual function to already be in terms of only xnew.
        #(i.e. plug in the value of the estimated derivative and time as well as the parameters). make in terms of xnew
        xold = x[:,n]
        residualsfromx!(res,xnew) = residualfunction(res,(xnew - xold)/dt,xnew,parameters,time)#x[n] grabs value of x for the current step.
        #Use the residual function with NLsolve to get xdot. use xold as starting guess.
        results=nlsolve(residualsfromx!,xold)#Solves equation and returns a "SolverResults" object.
        #store the result
        x[:,n+1] = results.zero #Grabs the zero (solution) field from the "SolverResults" object.
    end
    return x
end
#-----------End DISCRETE EULER
