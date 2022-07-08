"""
Эта программа соответствует коду, предложенному  и реализованному на языке С в статье
A Simple Numerical Algorithm and Software for Solution of Nucleation, Surface Growth, and Coagulation Problems
A. Prakash , A. P. Bapat & M. R. Zachariah (2003)

Суть: узловой метод решение общего уравнения динамики распределения частиц за счет процессов нуклеации, коагуляции и
гигроскопического роста/уменьшения размеров частиц (eq. 1).

Узловой метод менее затратный, чем метод деления диапазона на непересекающиеся интервалы.
В работе проведено тестирование модели в различных режимах, которое показало хорошее соответствие с имеющимися данными.

Athor: J. Kuznetsova
Date : 07.07.2022
"""
# загрузка функции-таймера
from timeit import default_timer as timer

# подключение библиотеки numpy
import numpy as np


start = timer()  # запуск таймера
end = timer()  # остановка таймера
print('время выполнения: %.3e с' % (end - start))
