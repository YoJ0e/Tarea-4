# ode/ode.py

""" # ODE 
Este módulo contiene varios métodos númericos para resolver ecuaciones
diferenciales ordinarias (ODEs por sus siglas en inglés). Permite al usuario
resolver la ODE que desee por medio de distintas funciones que implementan diversos
métodos.

El módulo contiene las siguientes funciones:

    - ´Euler´: Utiliza el Método de Euler para resolver la ODE.
    - ´RK2´: Utiliza el Método de RK2 para resolver la ODE.
    - ´RK4´: Utiliza el Método de RK4 para resolver la ODE.
"""

def Euler(f,xi,ti,tf,N):
    """Esta función utiliza el Método de Euler para resolver la ODE que el usuario dicte.

    Examples:
        >>> import numpy as np
        >>> def f(t,x):
        >>>     return -x**3-np.sin(t)
        >>> t,x = Euler(f,0,0,10,20)
        >>> print(t)
        [ 0.          0.52631579  1.05263158  1.57894737  2.10526316  2.63157895
        3.15789474  3.68421053  4.21052632  4.73684211  5.26315789  5.78947368
        6.31578947  6.84210526  7.36842105  7.89473684  8.42105263  8.94736842
        9.47368421 10.        ]
        >>> print(x)
        [ 0.          0.          0.26439534  0.71189381  1.04830651  0.89488942
        0.77464602  0.52141013  0.17502349 -0.28921312 -0.80263945 -0.97897518
        -0.73458315 -0.50879969 -0.16038526  0.30726691  0.81787725  0.97386719
	0.72957635  0.49945695]
		
    Args:
        f (función): Función de la ecuación que el usuario desee responder.
	xi (float | entero): Valor inicial para la variable dependiente.
        ti (float | entero): Valor inicial para la variable independiente.
        tf (float | entero): Valor final para la variable independiente.
        N (entero): Número de particiones entre el rango de variables independientes.

    Returns:
        t (vector): Un vector que contiene los resultados para las variables independientes 't'.
        x (vector): Un vector que contiene los resultados para las variables dependientes 'x'.
    """
    t = np.linspace(ti,tf,N)
    x = np.zeros(t.size)
    x[0] = xi
    h = t[1] - t[0]
    for i in range(t.size - 1):
        x[i+1] = x[i] + h*f(x[i],t[i])
    return t,x

def RK2(f,xi,ti,tf,N):
    """Esta función utiliza el Método de RK2 para resolver la ODE que el usuario requiera.

    Examples:
        >>> import numpy as np
        >>> def f(t,x):
        >>>     return -x**3-np.sin(t)
        >>> t,x = RK2(f,0,0,10,20)
        >>> print(t)
        [ 0.          0.52631579  1.05263158  1.57894737  2.10526316  2.63157895
        3.15789474  3.68421053  4.21052632  4.73684211  5.26315789  5.78947368
        6.31578947  6.84210526  7.36842105  7.89473684  8.42105263  8.94736842
        9.47368421 10.        ]
        >>> print(x)
        [ 0.          0.13691106  0.50040599  0.83221858  0.89696773  0.83638552
        0.684369    0.42791827  0.0377284  -0.46988042 -0.7896376  -0.78706438
        -0.65422559 -0.40234264 -0.0089812   0.49847514  0.79691697  0.78634182
        0.64914835  0.39298319]

    Args:
        f (función): Función de la ecuación que el usuario desee responder.
        xi (float | entero): Valor inicial para la variable dependiente.
        ti (float | entero): Valor inicial para la variable independiente.
        tf (float | entero): Valor final para la variable independiente.
        N (entero): Número de particiones entre el rango de variables independientes.

    Returns:
        t (vector): Un vector que contiene los resultados para las variables independientes 't'.
        x (vector): Un vector que contiene los resultados para las variables dependientes 'x'.
    """
    t = np.linspace(ti,tf,N)
    x = np.zeros(t.size)
    x[0] = xi
    h = t[1] - t[0]
    for i in range(t.size-1):
        k1 = h*f(x[i],t[i])
        k2 = h*f(x[i] + (k1/2), t[i] + (h/2))
        n = x[i] + (h/2)*f(x[i],t[i])
        x[i+1] = x[i] + h*f(n, t[i] + h/2)
    return t,x

def RK4(f,xi,ti,tf,N):
    """Esta función utiliza el Método de RK4 para resolver la ODE que el usuario brinde.

    Examples:
        >>> import numpy as np
        >>> def f(t,x):
        >>>     return -x**3-np.sin(t)
        >>> t,x = RK4(f,0,0,10,20)
        >>> print(t)
        [ 0.          0.52631579  1.05263158  1.57894737  2.10526316  2.63157895
        3.15789474  3.68421053  4.21052632  4.73684211  5.26315789  5.78947368
        6.31578947  6.84210526  7.36842105  7.89473684  8.42105263  8.94736842
        9.47368421 10.        ]
        >>> print(x)
        [ 0.          0.13505937  0.48487586  0.81979626  0.93102462  0.88144889
        0.72370945  0.46026479  0.06848691 -0.42773679 -0.78084333 -0.83096351
        -0.6993985  -0.43991777 -0.04506765  0.45089026  0.79057875  0.83047637
        0.69376603  0.43014634]

    Args:
        f (función): Función de la ecuación que el usuario desee responder.
        xi (float | entero): Valor inicial para la variable dependiente.
        ti (float | entero): Valor inicial para la variable independiente.
        tf (float | entero): Valor final para la variable independiente.
        N (entero): Número de particiones entre el rango de variables independientes.

    Returns:
        t (vector): Un vector que contiene los resultados para las variables independientes 't'.
        x (vector): Un vector que contiene los resultados para las variables dependientes 'x'.
    """
    t = np.linspace(ti,tf,N)
    x = np.zeros(t.size)
    x[0] = xi
    h = t[1] - t[0]
    for i in range(t.size-1):
        k1 = h*f(x[i],t[i])
        k2 = h*f(x[i] + (k1/2), t[i] + (h/2))
        k3 = h*f(x[i] + (k2/2), t[i] + (h/2))
        k4 = h*f(x[i] + k3, t[i] + h)
        x[i+1] = x[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    return t,x

