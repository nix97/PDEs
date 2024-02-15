'''
Wave Equation Solver by Lukas Setiawan :
e-mail: lukassetiawan@yahoo.com
My works on https://bitbucket.org/nixz97/nix/downloads/
         or https://github.com/nix97


Wave Equation is Partial Differential Equation(PDE) of hyperbolic type.
The solution using Finite-Difference method.

To approximate the solution of Wave Equation Utt(x,t)=c^2 Uxx(x,t)
over R = {(x,t): 0 ≤ x ≤ a, 0 ≤ t ≤ b with U(0,t)=0, U(a,t)=0 for 0 ≤ t ≤ b
and U(x,0)=F(x), Ut(x,0)=G(x) for 0 ≤ x ≤ a .
m,n are dimensions of the grid,(for i=1 to n; for j=1 to m)
and c is wave equation constant.
'''

from pymep.realParser import eval
import numpy as np
import plotly.graph_objects as go

again = "Y"
while again[0] in ("y", "Y"):

    print("Wave Equation Solver using Finite-Difference method in Python")
    print("Wave Equation Utt(x,t)=c^2 Uxx(x,t)")
    print("0 ≤ x ≤ a")
    print("0 ≤ t ≤ b")
    print("For example : Utt (x,t) = 4 Uxx (x,t); for 0 ≤ x ≤ 1  and  0 ≤ t ≤ 0.5")
    print("n=51 (dimension x); m=51 (dimension t)")
    print("F(x)=sin(pi*x)+sin(2*pi*x)")
    print("G(x)=0")
    print("From data and equation above we got a=1; b=0.5; c=2; n=51; m=51")
    print("F(x)=sin(pi*x)+sin(2*pi*x)")
    print("G(x)=0")
    print("a=1; b=0.5; c=2; m=51; n=51\n")

    print("Input:")
    F=input("F(x)=")
    G=input("G(x)=")
    a=float(input("a="))
    b=float(input("b="))
    c=float(input("c="))
    m=n=int(input("m=n="))

    h=a/(n-1)
    k=b/(m-1)
    r=(c*k)/h
    r2=r*r
    r22=(r*r)/2
    s1=1-(r*r)
    s2=2-(2*r*r)

    F3 = np.zeros(m+1)
    G3 = np.zeros(n+1)
    U = np.zeros((n+1,m+1))
    Z2=np.zeros((n+1,m+1))

    for i in range(2,n):
     x=h*(i-1)
     var = {"x": x}
     F3[i] =eval(F, var)
     G3[i]=eval(G, var)

    for j in range(1, m+1):
       U[1,j] = 0
       U[m,j] = 0


    for i in range(2,n):
     U[i,1]=F3[i]
     U[i,2] = s1*F3[i]+k*G3[i]+(r22* (F3[i + 1]+ F3[(i - 1)]))

    for j in range(3, m+1):
         for i in range(2, n ):
             U[i,j] = s2 * U[i,j - 1] + r2 * (U[i - 1,j - 1] + U[i + 1,j - 1]) - U[i,j - 2]

    print("")
    print("Result:")
    print("{:<8} {:<10} {:<10} {:<10} {:<30}".format('i', 'j', 'Xi', 'Tj','U[i,j]'))

    for j in range(1,m+1):
     for i in range(1,n+1):
       x = h * (i - 1)
       t = k * (j - 1)
       print ("{:<8} {:<10} {:<10.3} {:<10.3} {:<30}".format(i,j,x, t,U[i,j]))

    sh_0, sh_1 = n,m
    x, t = np.linspace(0, a, n), np.linspace(0, b, m)

    #Display Graph 3D(surface) on browser
    fig = go.Figure(data=[go.Surface(z=U, x=x, y=t)])
    fig.update_layout(
        title='Graphic Wave Equation',
        autosize=True,
        scene=dict(
            xaxis_title='X',
            yaxis_title='T',
            zaxis_title='U',
        ),
    )

    print("Wait the moment to show graph on browser...")
    fig.show()

    again = input("Another one (Y/N)?")

