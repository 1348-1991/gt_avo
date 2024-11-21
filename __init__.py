# -*- coding: utf-8 -*-
#    ===================
#    gt_avo.py
#    ===================
# Geophysical Tools: AVO
# Conjunto de Funções coletadas para Análise de AVO e afins.
#
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

from . import gt_avo
from . import akirichards
from . import shuey
from . import shuey2
from . import intercept
from . import gradient
from . import earth_model
from . import ricker_wav
from . import ref_c
from . import a_gather
from . import off_gather
from . import norm_wavelet
