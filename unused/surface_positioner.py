# -*- coding: utf-8 -*-
######################################################################
# This file is part of sage-flatsurf.
#
#       Copyright (C) 2016 Pat Hooper
#                     2020 Julian Rüth
#
# sage-flatsurf is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# sage-flatsurf is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with sage-flatsurf. If not, see <https://www.gnu.org/licenses/>.
######################################################################
from sage.all import Fields

class SurfacePositioner:
    def __init__(self,field):
        if not field in Fields():
            raise TypeError("field must be a field")

        self._positions={}
        self._field=field
        self._zero=self.vector_space().zero()

    def get_position(self,index):
        try:
            pos=self._positions[index]
        except KeyError:
            pos=self._zero
        return pos

    def translate_position(self, index, translation):
        curpos=self.getPosition(index)
        try:
            newpos=curpos+translation
        except StandardError:
            from sage.modules.free_module_element import vector
            try:
                newpos=curpos + vector(translation)
            except StandardError:
                try:
                    newpos=curpos+vector(self._field(translation[0]),self._field(translation[1]))
                except:
                    raise ValueError("can not convert translation to a vector in R^2")
        if newpos==self._zero:
            try:
                del self._positions[index]
            except KeyError:
                pass
        else:
            self._positions[index]=newpos

    def set_position(self, index, position):
        try:
            position=self._zero+position
        except StandardError:
            from sage.modules.free_module_element import vector
            try:
                position=self._zero+vector(position)
            except StandardError:
                try:
                    position=vector(self._field(translation[0]),self._field(translation[1]))
                except:
                    raise ValueError("can not convert translation to a vector in R^2")
        if position==self._zero:
            try:
                del self._positions[index]
            except KeyError:
                pass
        else:
            self._positions[index]=position

    def vector_space(self):
        r"""
        Return the vector space in which self naturally embeds.
        """
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self._field, 2)

