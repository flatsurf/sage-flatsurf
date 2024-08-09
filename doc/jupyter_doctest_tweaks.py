r"""
Sanitizes jupyter-sphinx cells in the documentation.

We use jupyter-execute to render plots in the module reference manual. At the
same time, we want to use these jupyter-execute blocks as doctests. However,
the sample output in these doctests cause syntax errors when passed to the
SageMath Jupyter kernel.

Here, we drop all lines that do not start with sage: or ....: so that the
Jupyter kernel can execute them.
"""
# ****************************************************************************
#  This file is part of sage-flatsurf.
#
#       Copyright (C) 2024 Julian Rüth <julian.rueth@fsfe.org>
#
#  sage-flatsurf is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  sage-flatsurf is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sage-flatsurf. If not, see <https://www.gnu.org/licenses/>.
# ****************************************************************************

import re
from IPython.core.inputtransformer2 import PromptStripper
import IPython

# Replace the transform that strips sage: and ...: with our more aggressive version.
IPython.get_ipython().input_transformer_manager.cleanup_transforms[1] = PromptStripper(
    prompt_re=re.compile(r"^((sage: |\.\.\.\.: )|(?!sage: |\.\.\.\.: ).*)")
)
