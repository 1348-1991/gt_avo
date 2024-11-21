"""
===================
gt_avo.py
===================
Geophysical Tools: AVO
Conjunto de Funções coletadas para Análise de AVO e afins.

"""

# A Python Module for Manipulating AVO related operations.
# Created November 2024 by Vitor Azevedo dos Santos (vitors@id.uff.br)
# History: 2024-11-21 First Public Release.
# Copyright: n/a
# License: n/a
# ++ updates

# Main Functions:

# gt_avo.akirichards                      : Aproximação de Aki Richards
# gt_avo.shuey                            : Aproximação de Shuey com 3 Termos
# gt_avo.shuey2                           : Aproximação de Shuey com 2 Termos
# gt_avo.intercept                        : Cálculo de Intercept
# gt_avo.gradient                         : Cálculo de Gradient
# gt_avo.earth_model                      : Modelagem da Terra
# gt_avo.ricker_wav                       : Wavelet de Ricker
# gt_avo.ref_c                            : Coeficiente de Reflexão
# gt_avo.a_gather                         : Angle Gather
# gt_avo.off_gather                       : Offset Gather
# gt_avo.norm_wavelet                     : Normalização de Wavelet

#   THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
#   WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
#   MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
#   ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
#   OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
#   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
#   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
#   ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Definindo Aki-Richard (função: akirichards):

def akirichards(vp1, vs1, rho1, vp2, vs2, rho2, theta1):

  """
  Calcula a resposta AVO para reflexão baseada em Aki-Richard em 3 termos,
  aproximação das equações de Zoeppritz.

  Args:
    vp1: Velocidade das ondas P do meio superior (m/s).
    vs1: Velocidade das ondas S do meio superior (m/s).
    rho1: Densidade do meio superior (kg/m^3).

    vp2: Velocidade das ondas P do meio inferior (m/s).
    vs2: Velocidade das ondas S do meio inferior (m/s).
    rho2: Densidade do meio inferior (kg/m^3).
    theta1: Ângulo de incidência da onda P na camada superior (em graus).

    Retorna:
    R: Coeficiente de Reflexão.
  """
  # Converte o ângulo de incidência para radianos
  theta_rad = np.deg2rad(theta1)

  # Calcula os termos intermediários:
  # Calcula dos termos intermediários: variações e médias
  dvp = vp2 - vp1                 # delta vp
  dvs = vs2 - vs1                 # delta vs
  drho = rho2 - rho1              # delta rho
  mvp = (vp1 + vp2) / 2.0         # média vp
  mvs = (vs1 + vs2) / 2.0         # média vs
  mrho = (rho1 + rho2) / 2.0      # média rho

  # Cálculo de Refletividade com 3 Termos:
  # Calcula o Coeficiente de Reflexão de Aki-Richards (1979): R
  # Calcula os termos em A, B, C  e p:
  p = np.sin(theta_rad) / vp1
  A = + 0.5 * (1 - 4 * p**2 * mvs**2) * (drho / mrho)
  B = + (dvp / (2 * np.cos(theta_rad)**2 * mvp))
  C = - (4 * p**2 * mvs**2 * dvs / mvs)

  R_Aki_Richard = A + B + C

  return R_Aki_Richard

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Definindo Shuey com 3 termos (função: shuey):

def shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
  """
  Calcula a resposta AVO para reflexão baseada em Shuey em 3 termos,
  aproximação das equações de Zoeppritz.

    Args:
    vp1: Velocidade das ondas P do meio superior (m/s).
    vs1: Velocidade das ondas S do meio superior (m/s).
    rho1: Densidade do meio superior (kg/m^3).

    vp2: Velocidade das ondas P do meio inferior (m/s).
    vs2: Velocidade das ondas S do meio inferior (m/s).
    rho2: Densidade do meio inferior (kg/m^3).

    theta1: Ângulo de incidência da onda P na camada superior (em graus).

    Retorna:
    R: Coeficiente de Reflexão.
  """
  # Converte o ângulo de incidência para radianos:
  theta_rad = np.deg2rad(theta1)

  # Calcula dos termos intermediários: variações e médias
  dvp = vp2 - vp1                 # delta vp
  dvs = vs2 - vs1                 # delta vs
  drho = rho2 - rho1              # delta rho
  mvp = (vp1 + vp2) / 2.0         # média vp
  mvs = (vs1 + vs2) / 2.0         # média vs
  mrho = (rho1 + rho2) / 2.0      # média rho

  # Cálculo de Refletividade com 3 Termos:
  # Calcula o Coeficiente de Reflexão Shuey a partir de Aki-Richards (Conforme Shuey 1985): R
  # Calcula os termos em A, B, and C:

  A = 0.5 * (dvp / mvp + drho / mrho) # Interceptação (ângulo normal)
  B = 0.5 * dvp / mvp - 2 * mvs**2 / mvp**2 * (drho / mrho + 2 * dvs / mvs)
  C = 0.5 * dvp / mvp

  R_Shuey_3 = A + B * np.sin(theta_rad)**2 + C * (np.tan(theta_rad)**2 - np.sin(theta_rad)**2)

  return R_Shuey_3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Definindo Shuey com 2 termos (função: shuey2):

def shuey2(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
  """
  Calcula a resposta AVO para reflexão baseada em Shuey em 2 termos,
  aproximação das equações de Zoeppritz.

    Args:
    vp1: Velocidade das ondas P do meio superior (m/s).
    vs1: Velocidade das ondas S do meio superior (m/s).
    rho1: Densidade do meio superior (kg/m^3).

    vp2: Velocidade das ondas P do meio inferior (m/s).
    vs2: Velocidade das ondas S do meio inferior (m/s).
    rho2: Densidade do meio inferior (kg/m^3).

    theta1: Ângulo de incidência da onda P na camada superior (em graus).

    Retorna:
    R: Coeficiente de Reflexão.
  """
  # Converte o ângulo de incidência para radianos:
  theta_rad = np.deg2rad(theta1)

  # Calcula dos termos intermediários: variações e médias
  dvp = vp2 - vp1                 # delta vp
  dvs = vs2 - vs1                 # delta vs
  drho = rho2 - rho1              # delta rho
  mvp = (vp1 + vp2) / 2.0         # média vp
  mvs = (vs1 + vs2) / 2.0         # média vs
  mrho = (rho1 + rho2) / 2.0      # média rho

  # Cálculo de Refletividade com 2 Termos:
  # Calcula o Coeficiente de Reflexão Shuey a partir de Aki-Richards (Conforme Shuey 1985): R
  # Calcula os termos em A e B:

  A = 0.5 * (dvp / mvp + drho / mrho) # Interceptação (ângulo normal)
  B = ((0.5*(dvp / mvp) - 4*(mvs**2 / mvp**2) * (dvs / mvs) - 2*(mvs**2 / mvp**2) * (drho / mrho)) * np.sin(theta_rad)**2) # Ângulos Intermediários

  R_Shuey_2 = A + B

  return R_Shuey_2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Definindo Intercept (função: intercept)

def intercept(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
  """
  Calcula o Intercept baseado em Aki-Richard,
  aproximação das equações de Zoeppritz.

  Args:
    vp1: Velocidade das ondas P do meio superior (m/s).
    vs1: Velocidade das ondas S do meio superior (m/s).
    rho1: Densidade do meio superior (kg/m^3).

    vp2: Velocidade das ondas P do meio inferior (m/s).
    vs2: Velocidade das ondas S do meio inferior (m/s).
    rho2: Densidade do meio inferior (kg/m^3).
    theta1: Ângulo de incidência (graus).

    Retorna:
    I: Intercept
  """
  # Converte o ângulo de incidência para radianos
  theta_rad = np.deg2rad(theta1)

  # Calcula os termos intermediários:
  # Calcula dos termos intermediários: variações e médias
  dvp = vp2 - vp1                 # delta vp
  dvs = vs2 - vs1                 # delta vs
  drho = rho2 - rho1              # delta rho
  mvp = (vp1 + vp2) / 2.0         # média vp
  mvs = (vs1 + vs2) / 2.0         # média vs
  mrho = (rho1 + rho2) / 2.0      # média rho

  # Calcula o Intercept: I

  #intercepto = 0.5 * (dvp / mvp + drho / mrho) # Interceptação (ângulo normal)
  #G = 2 * (mvs**2 / mvp**2) * (2 * (dvs / mvs + drho / mrho))
  #intercepto = 0.5 * ((dvp / mvp) + (drho / mrho) + G * np.sin(theta_rad)**2)

  intercepto = 0.5 * ((dvp / mvp) + (drho / mrho))

  return intercepto

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Definindo Gradient (função gradient):

def gradient(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
  """
  Calcula o Gradient
  Args:
    vp1: Velocidade das ondas P do meio superior (m/s).
    vs1: Velocidade das ondas S do meio superior (m/s).
    rho1: Densidade do meio superior (kg/m^3).

    vp2: Velocidade das ondas P do meio inferior (m/s).
    vs2: Velocidade das ondas S do meio inferior (m/s).
    rho2: Densidade do meio inferior (kg/m^3).
    theta1: Ângulo de incidência (graus).

    Retorna:
    G: Gradiente
  """
  # Converte o ângulo de incidência para radianos
  theta_rad = np.deg2rad(theta1)

  # Calcula os termos intermediários:
  # Calcula dos termos intermediários: variações e médias
  dvp = vp2 - vp1                 # delta vp
  dvs = vs2 - vs1                 # delta vs
  drho = rho2 - rho1              # delta rho
  mvp = (vp1 + vp2) / 2.0         # média vp
  mvs = (vs1 + vs2) / 2.0         # média vs
  mrho = (rho1 + rho2) / 2.0      # média rho

  # Calcula o Gradient: G
  # G_Gradiente = (2*(mvs**2 / mvp**2) * (2*(dvs / mvs + drho / mrho))) * np.sin(theta_rad)**2 # Gradiente
  G_Gradiente = 2 * (mvs**2 / mvp**2) * (2 * (dvs / mvs + drho / mrho))

  return G_Gradiente

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Criando Modelo da terra 1 D (função earth_model):

def earth_model(vp1, vs1, rho1, vp2, vs2, rho2, n_samples):
    """Cria um modelo de Terra 1D com duas camadas.

    Args:
        vp1, vs1, rho1: Propriedades da primeira camada.
        vp2, vs2, rho2: Propriedades da segunda camada.
        n_samples: Número de amostras.

    Returns:
        np.ndarray: Modelo de velocidade P.
        np.ndarray: Modelo de velocidade S.
        np.ndarray: Modelo de densidade.
    """

    vp = np.concatenate(([vp1] * int(n_samples/3), [vp2] * int(n_samples/3), [vp1] * int(n_samples/3)))
    vs = np.concatenate(([vs1] * int(n_samples/3), [vs2] * int(n_samples/3), [vs1] * int(n_samples/3)))
    rho = np.concatenate(([rho1] * int(n_samples/3), [rho2] * int(n_samples/3), [rho1] * int(n_samples/3)))

    return vp, vs, rho

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Definindo Wavelet de Ricker (função ricker_wav):

def ricker_wav(dt, freq, length):
  """
  Gera uma wavelet de Ricker.

  Args:
    dt: Intervalo de tempo entre amostras.
    freq: Frequência dominante da wavelet.
    length: Duração total da wavelet em segundos.

  Returns:
    time: Array com os valores de tempo.
    wavelet: Array com os valores da wavelet.
  """

  t = np.arange(0, length, dt)
  sigma = 1 / (np.pi * freq)
  wavelet = (1 - 2 * (np.pi**2) * (freq**2) * (t**2)) * np.exp(-(np.pi**2) * (freq**2) * (t**2))

  # Normalização
  # wavelet /= np.max(np.abs(wavelet))

  return t, wavelet # Retorna Tempo e a Wavelet

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Definindo Coeficiente de Reflexão (função ref_c):

def ref_c(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
    """Calcula os coeficientes de reflexão para diferentes ângulos de incidência.

    Args:
        vp1, vs1, rho1: Propriedades da primeira camada.
        vp2, vs2, rho2: Propriedades da segunda camada.
        angles: Ângulos de incidência em radianos.

    Returns:
        tuple: (reflect, r0, g) onde:
            reflect: Coeficientes de reflexão para os ângulos fornecidos.
            r0: Coeficiente de reflexão no zero-offset.
            g: Gradiente.
    """
    # Conversão para radianos:
    theta_rad = np.deg2rad(theta1)

    # Impedância Acústica das Camadas
    ip1 = vp1 * rho1
    ip2 = vp2 * rho2

    # Coeficientes de reflexão para Incidência Normal (R0) e Gradiente (G)
    R0 = (ip2 - ip1) / (ip2 + ip1)
    G = 0.5 * (1 - (4 * vs1**2 / vp1**2)) * (1 - (vp2 / vp1)**2)

    # Aproximação de Shuey
    refco = R0 + G * np.sin(theta_rad)**2

    # Retorna o Coeficiente de Reflexão (refco), o em Zero Offset (R0) e o Gradiente (G)
    return refco, R0, G

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Definindo Angle Gather (função a_gather):

def a_gather(rc, wavelet, angle_range):
    """
    Gera um Angle Gather sintético.

    Args:
        rc: Coeficientes de Reflexão, matriz 2D (Profundidade x Ângulos).
        wavelet: Wavelet Ricker para Convolução.
        angle_range: Faixa de ângulos.

    Returns:
        np.ndarray: Angle Gather no formato (Profundidade x Ângulos).
    """
    # Conversão para radianos:
    theta_rad = np.deg2rad(angle_range)

    # Lista para armazenar os resultados por ângulo
    anglegather = []

    # rc precisa ser um array 2D com dimensões (número de amostras, número de ângulos)
    if rc.ndim == 1:
        # Se rc for 1D, assumimos que contém coeficientes para um único ângulo
        # e replicamos para todos os ângulos no angle_range
        num_samples = len(rc)
        num_angles = len(angle_range)
        rc = np.tile(rc[:, np.newaxis], (1, num_angles))  # Replica rc para todos os ângulos

    for i in range(len(angle_range)):
        # Convolução da wavelet com os coeficientes de reflexão para o ângulo i
        anglegather.append(np.convolve(wavelet, rc[:, i], mode='same'))

    # Converte para array e transpõe para (Tempo (Profundidade) , Ângulos)
    return np.array(anglegather).T

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Definindo Offset Gather (função off_gather):

def off_gather(rc, wavelet, offset_range, angle_range, velocity):
    """Cria um offset gather.

    Args:
        rc: Coeficientes de Reflexão.
        wavelet: Wavelet Ricker.
        offset_range: Faixa de offsets.
        angle_range: Faixa de ângulos de incidência.
        velocity: Velocidade da onda P na primeira camada (m/s).
        velocity: Importante para o cálculo do Delay

    Returns:
        np.ndarray: Offset gather.
    """
    # Conversão para radianos:
    angle_range_rad = np.deg2rad(angle_range)

    gather = []

    # Verifica as dimensões de rc e replica para cada ângulo se necessário:
    if rc.ndim == 1:
        num_angles = len(angle_range_rad)
        rc = np.tile(rc[:, np.newaxis], (1, num_angles))

    # Loop para cada Offset e Ângulo:
    for offset in offset_range:
        for angle_index in range(len(angle_range_rad)):
            # Simulação de Atraso usando a velocidade (em geral vp1):
            delay = int(offset / velocity)

            # Convolução com atraso, usando rc para o ângulo atual:
            trace = np.convolve(rc[:, angle_index], wavelet, mode='same')
            trace = np.roll(trace, delay)
            gather.append(trace)

    # Reorganiza o gather para as dimensões (tempo, offset, ângulo)
    gather = np.array(gather).reshape(len(offset_range), len(angle_range_rad), -1)
    gather = gather.transpose(2, 0, 1)  # Transpõe para (Tempo, Offset, Ângulo)

    return gather

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Normalizando Wavelet de Ricker (função norm_wavelet):

def norm_wavelet(wavelet):
    """Normaliza uma wavelet para que sua amplitude máxima seja 1.

    Args:
        wavelet: A wavelet a ser normalizada.

    Returns:
        A wavelet normalizada.
    """

    amplitude_maxima = np.max(np.abs(wavelet))
    norm_factor = 1 / amplitude_maxima
    wavelet_norm = wavelet * norm_factor
    return wavelet_norm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
