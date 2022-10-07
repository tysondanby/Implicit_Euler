
function DAEroberts!(res,xdot,x,params,time)
    res[1] = - 0.04x[1]              + 1e4*x[2]*x[3] - xdot[1]
    res[2] = + 0.04x[1] - 3e7*x[2]^2 - 1e4*x[2]*x[3] - xdot[2]
    res[3] = x[1] + x[2] + x[3] - 1.0

end

function DAElinear!(res,xdot,x,params,time)
    res[1] = 2*xdot[1] + .8*xdot[2]
    res[2] = -1*xdot[1] + -3*xdot[2]
    res[3] = 5-xdot[3]
end

function DAEpendulum!(res,xdot,x,params,time)
    th1,th2,x2,y2,om1,om2,u,v,T = x
    dth1,dth2,dx2,dy2,dom1,dom2,du,dv,dT = x
    m1,m2,l1,l2,g = params

    #torque balance, top bob
    res[1]= -dom1*m1*l1 -sin(th1)*m1*g + T*sin(th2)
    #force balance, bottom bob
    res[2]= du*m2 +sin(th1+th2)*T
    res[3]= -dv*m2 +cos(th1+th2)*T - m2*g
    #Geometric constraint for relating th2 to th1, x, and y.
    #res[9]= atan(x2-l1*sin(th1),-y2-l1*cos(th1)) - th1 - th2
    #Geometric constraint for making sure the length of the second bob is l2
    res[4]= l2^2 - (x2-l1*sin(th1))^2 - (-y2-l1*cos(th1))^2
    #Some state variables are the derivatives of others.
    res[5:8]=xdot[1:4]-x[5:8]
end

function simplependulum!(res,xdot,x,params,time)
    th1,th2,om1,om2 = x
    dth1,dth2,dom1,dom2 = xdot
    m1,m2,l1,l2,g = params

    res[1] = -dom1 + (-g*(2*m1 + m2)*sin(th1) - m2*g*sin(th1 - 2*th2) - 2*sin(th1 - th2)*m2*(om2^2*l2 + om1^2*l1*cos(th1 - th2)))/(l1*(2*m1 + m2 - m2*cos(2*th1 - 2*th2)))
    res[2] = -dom2 + (2*sin(th1-th2)*(om1^2*l1*(m1 + m2) + g*(m1 + m2)*cos(th1) + om2^2*l2*m2*cos(th1 - th2)))/(l2*(2*m1 + m2 - m2*cos(2*th1 - 2*th2)))
    res[3] = om1 - dth1
    res[4] = om2 - dth2
end
