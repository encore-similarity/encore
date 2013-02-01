# utils.py --- Utility functions for ensemble package
# Copyright (C) 2012 Wouter Boomsma
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


def vararg_callback(option, opt_str, value, parser):
    '''A callback for the option parser allowing a variable number of arguments.'''

    value = []

    for arg in parser.rargs:

        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break

        # stop on -a, but not on -3 or -3.0
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break

        value.append(arg)
    del parser.rargs[:len(value)]

    # parser.values.ensure_value(option.dest, []).append(value)
    setattr(parser.values, option.dest, value)
