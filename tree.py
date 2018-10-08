# Tree Generation Tool
# Copyright (C) 2018 Alexander Kraus <nr4@z10.info>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
 

import numpy as np
import matplotlib.pyplot as plt
import sys

initiator = "X"
state = initiator
depth = 5

for i in range(depth):
    state = state.replace('X', 'F+[[X]-X]-F[-FX]+X')
    state = state.replace('F', 'FF')
    print(state)

tpos = [ 0., 0., .5*np.pi]
stack = []
lines = []
xmin = sys.maxsize
xmax = -sys.maxsize
ymin = sys.maxsize
ymax = -sys.maxsize

chains = state.replace('[',',').replace(']',',').replace('+',',').replace('-',',').replace('X',',').split(',')
chains = list(set(chains))
lens = [ len(chain) for chain in chains ]

simple = state
for leni in reversed(sorted(lens)):
    print(leni)
    simple = simple.replace('F'*leni, chr(ord('a')+int(np.log2(leni+1))))
print("simple:", simple)

print(sorted(lens))

for c in state:
    xmin = min(xmin, tpos[0])
    xmax = max(xmax, tpos[0])
    ymin = min(ymin, tpos[1])
    ymax = max(ymax, tpos[1])
        
    if c == 'F': # Line forward
        p0 = [ tpos[0], tpos[1] ]
        tpos = [tpos[0] + np.cos(tpos[2]), tpos[1] + np.sin(tpos[2]), tpos[2]]
        p1 = [ tpos[0], tpos[1] ]
        lines += [ [ p0, p1 ] ]
    elif c == ']':
        tpos = stack.pop()
    elif c == '[':
        stack += [ tpos ]
    elif c == '+':
        tpos = [tpos[0], tpos[1], tpos[2]+.001*np.pi]
    elif c == '-':
        tpos = [tpos[0], tpos[1], tpos[2]-.001*np.pi]

for i in range(len(lines)):
    for j in range(len(lines[i])):
        lines[i][j][0] /= abs(xmax-xmin)
        lines[i][j][1] /= abs(ymax-ymin)

fig = plt.figure()
for line in lines:
    p0 = line[0]
    p1 = line[1]
    
    x = [ p0[0]*(1.-t)+p1[0]*t for t in np.arange(0., 1.1, 1./3.) ]
    y = [ p0[1]*(1.-t)+p1[1]*t for t in np.arange(0., 1.1, 1./3.) ]
    
    plt.plot(x, y, '-')
    
plt.show()
