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
    return x,t

def RK4(f,xi,ti,tf,N):
    """Esta función utiliza el Método de RK4 para resolver la ODE que el usuario brinde.

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
    return x,t

