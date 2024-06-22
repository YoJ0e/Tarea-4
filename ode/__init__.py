# ode/__init__.py

"""
En este módulo se van a encontrar funciones que
describen varios métodos númericos para resolver ecuaciones
diferenciales ordinarias (ODEs por sus siglas en inglés).

El módulo exportado con este paquete es:

    - ´ode´: Contiene dentro del archivo 'ode.py' los métodos númericos de Euler, 
             RK2 y RK4. Estos pueden ser implementados en diferentes funciones 
             por medio de su argumento para resolver la ODE.

## Método de Euler

Este método se basa en la expansión de Taylor de una función $x(t)$, el cual va a ser:
$$
x(t+h) = x(t) + h\\frac{dx}{dt} + \\overbrace{ \\frac{h^2}{2} \\frac{d^2x}{dt^2} } ^{\\epsilon} + O(h^3)
$$
donde $\\epsilon$ es el error asociado al cálculo, el cual no se va a calcular. Este error depende de la cantidad
de veces que se divide el intervalo de variables independientes.

De aquí se puede ver que para avanzar una cantidad $\\textbf{h}$ entre el intervalo
de variables independientes basta con solo utilizar la ecuación:
$$
x(t + h) = x(t) + hf(x,t)
$$
donde $t$ es la variable independiente y $x$ la dependiente.

Cabe mencionar que para realizar el cálculo, se tiene que dividir el intervalo de variables independientes
en N secciones, por lo que hace falta hacer el cálculo $N = (b-a)/h$. Con todo eso en cuenta se encuentra que
la ecuación necesaria para utilizar el método de Euler va a ser:
$$
x_i = x_{i-1} + hf(x_{i-1})
$$

## Método de RK2

El método de Runge-Kutta es una familia de métodos de distinto orden que se pueden utilizar para mejorar la
aproximación realizada en el método de Euler sin necesidad de considerar ordenes mayores en la expansión de
Taylor. En este módulo inicialmente se implementó el método de Runge-Kutta de segundo orden (RK2). Este método
consiste en realizar un cálculo muy similar al de Euler, pero utilizando el punto medio entre los saltos para
realizar la aproximación, mejorando así la aproximación que se va a realizar con la misma $\\textbf{h}$. Tomando
en cuenta estas consideraciones, se llega a la ecuación:
$$
x(t + h) = x(t) + hf[x(t + h/2), t + h/2] + O(h^3)
$$

Sin embargo surge un problema, el hecho de que no se puede calcular el valor de la función en el punto medio
$x(t+h/2)$. Para aproximar este valor, se utiliza el método de Euler con un paso $h/2$. De esta manera se encuentran las ecuaciones:
$$
k_1 = hf(x,t)
$$
$$
k_2 = hf\\left(x + \\frac{k_1}{2},t + \\frac{h}{2}\\right)
$$
$$
x(t+h) = x(t) + k_2
$$

Cabe mencionar que este cálculo posee el mismo error que el realizado con el Método de Euler.

## Método de RK4

Nuevamente, este es un método de Runge-Kutta, solo que esta vez de orden 4. Consiste en aplicar la misma metodología anterior
pero para más puntos ubicados entre $x(t)$ y $x(t+h)$. El problema con el Método de Runge-Kutta es que se encuentran
expresiones más complicadas conforme se va aumentando el orden de la aproximación. Es por esto que RK4 corresponde al mejor
balance entre complejidad y aproximación en esta familia de métodos. Este es el método más común a la hora de resolver ODEs.

Para aplicar este método se requiere utilizar las ecuaciones:
$$
k_1 = hf(x, t)
$$
$$
k_2 = hf\\left(x + \\frac{k_1}{2}, t+\\frac{h}2\\right)
$$
$$
k_3 = hf\\left(x + \\frac{k_2}{2}, t+\\frac{h}2\\right)
$$
$$
k_4 = hf\\left(x + k_3, t + h \\right)
$$
$$
x(t+h) = x(t) + \\frac{1}{6}(k_1 + 2 k_2 + 2k_3 + k_4)
$$

El error de aproximación asociado a este método es del orden de $O(h^5)$.
"""
