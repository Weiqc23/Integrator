"""PyWENO quadrature points.

Requires SymPy.
Namdi: I have edited this code from the original PyWENO code.
"""

# Copyright (c) 2011, Matthew Emmett.  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#   1. Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#   2. Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.



# ===============================================
# import
# ===============================================

import sympy
#from sympy import mpmath
import mpmath
################################################################################
# polynomial generator, roots etc

# ===============================================
# constants
# ===============================================
GAUSS_LEGENDRE  = 1
GAUSS_LOBATTO   = 2
GAUSS_RADAU     = 3
GAUSS_RADAU_2A  = 4
#UNIFORM = 5


# ===============================================
# function
# ===============================================


def find_roots(p):
  """
  Return set of roots of polynomial *p*.

  This uses the *nroots* method of the SymPy polynomial class to give
  rough roots, and subsequently refines these roots to arbitrary
  precision using mpmath.

  :param p: sympy polynomial

  :return: sorted *set* of roots.
  """

  # 声明 SymPy 符号变量 x，用于后续多项式求值
  x = sympy.var('x')
  
  # 创建空集合存储根（利用集合特性自动去重，避免重复根）
  roots = set()

  # 遍历 SymPy 多项式 p 的粗略根估计（p.nroots() 返回多项式的数值根近似）
  for x0 in p.nroots():
    # 1. 定义匿名函数：输入 z，返回多项式 p 在 x=z 处的取值（p.eval(x, z) 计算多项式在 z 处的值）
    # 2. 调用 mpmath.findroot，以 x0 为初始猜测值，高精度求解方程 p(z)=0 的根
    xi = mpmath.findroot(lambda z: p.eval(x, z), x0)
    
    # 将精细求解的根添加到集合（自动去重）
    roots.add(xi)

  # 将集合转换为排序后的列表并返回（确保根按数值大小升序排列）
  return sorted(roots)

def gauss_legendre(n):

  """
  This function returns Gauss-Legendre nodes

  Gauss-Legendre nodes are roots of :math:`P_n(x)`.

  :param int n: the number of nodes

  :return: return Gauss-Legendre nodes :math:`x \\in [-1, 1]`
  """

  p = legendre_poly(n)
  r = find_roots(p)
  return r


def gauss_lobatto(n):
    """
    This function returns Gauss-Lobatto nodes. Gauss-Lobatto nodes are roots of :math:`P'_{n-1}(x)`.

    :param int n: the number of nodes
    :return: Gauss-Lobatto nodes :math:`x \\in [-1, 1]`
    """
    x = sympy.var('x')                  # 步骤1：声明符号变量 x
    p = legendre_poly(n-1).diff(x)      # 步骤2：计算 (n-1) 次 Legendre 多项式的导数
    r = find_roots(p)                   # 步骤3：求解导数多项式的根（内部节点）
    r = [mpmath.mpf('-1.0'), mpmath.mpf('1.0')] + r  # 步骤4：添加端点 [-1, 1]
    return sorted(r)                    # 步骤5：排序并返回节点


def gauss_radau(n):

  """
  Return Gauss-Radau nodes (right hand left hand point). Gauss-Radau nodes are roots of \
  :math:`P_n(x) + P_{n-1}(x)`.

  :param int n: the number of nodes

  :return: Gauss-Radau nodes :math:`x \\in [-1, 1]`
  """

  p = legendre_poly(n) + legendre_poly(n-1)
  r = find_roots(p)

  return r

def gauss_radau_2a(n):

  """
  Return Gauss-Radau 2a nodes (left hand left hand point). Gauss-Radau 2a nodes are roots of \
  :math:`P_n(x) - P_{n-1}(x)`.

  :param int n: the number of nodes

  :return: Gauss-Radau 2a nodes :math:`x \\in [-1, 1]`
  """

  p = legendre_poly(n) - legendre_poly(n-1)
  r = find_roots(p)

  return r

def legendre_poly(n):

  """
  Return Legendre polynomial :math:`P_n(x)`.

  :param int n: the degree of the polynomial

  :return: Legendre polynomial
  """

  x = sympy.var('x')
  p = (1.0*x**2 - 1.0)**n

  top = p.diff(x, n)
  bot = 2**n * 1.0*sympy.factorial(n)

  return (top / bot).as_poly()