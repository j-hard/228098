
from pyplasm import *
pts0=[[-20,-20],[73,-20],[73,73],[-20,73]]
pts1=[[0,0],[53,0],[53,53],[0,53]]
pts2=[[6.2,6.2],[46.8,6.2],[46.8,46.8],[6.2,46.8]]
cells=[[1,2,3,4]]
pols=None
base=MKPOL([pts0,cells,pols])
floor1=T([1,2,3])([0,0,10])(PROD([MKPOL([pts1,cells,pols]),Q(-0.3)]))
floor2=T([1,2,3])([0,0,19])(PROD([MKPOL([pts1,cells,pols]),Q(-0.3)]))
floor3=T([1,2,3])([0,0,28])(PROD([MKPOL([pts1,cells,pols]),Q(-0.3)]))
floor4=T([1,2,3])([0,0,37])(PROD([MKPOL([pts1,cells,pols]),Q(-0.3)]))
floor5=T([1,2,3])([0,0,46])(PROD([MKPOL([pts1,cells,pols]),Q(-0.3)]))
floor6=T([1,2,3])([0,0,55])(PROD([MKPOL([pts1,cells,pols]),Q(-0.3)]))
copertura=T([1,2,3])([0,0,64])(MKPOL([pts2,cells,pols]))
intermedie=STRUCT([floor1,floor2,floor3,floor4,floor5,floor6])
plants=STRUCT([base,intermedie,copertura])
pts0=[[0,0],[53,0],[53,60],[0,60]]
pts1=[[6.2,0],[46.8,0],[46.8,64],[6.2,64]]

cells0=[[1,2,3,4]]
pols0=None
pts3=[[5.8,0],[0.8,0],[0.8,9.7],[5.8,9.7]]
pts4=[[5.8,0],[0.8,0],[0.8,8.7],[5.8,8.7]]

window0=MKPOL([pts3,cells0,pols0])
windows0=STRUCT([T([1])([(i*5.8)])(window0) for i in range(9)])
plane0=DIFFERENCE([MKPOL([pts0,cells0,pols0]), windows0])

windowl1=MKPOL([pts4,cells0,pols0])
windowsl1=STRUCT([T([1,2])([(i*5.8),10])(windowl1) for i in range(9)])
windows1=STRUCT([T([2])([(i*9)])(windowsl1) for i in range(5)])
plane=DIFFERENCE([plane0, windows1])
#plane=plane0

windows01=STRUCT([T([1])([5.7+(i*5.8)])(window0) for i in range(7)])


windowsl01=STRUCT([T([1,2])([5.7+(i*5.8),10])(windowl1) for i in range(7)])
windows01=STRUCT([T([2])([(i*9)])(windowsl01) for i in range(5)])
windows00=STRUCT([T([1])([(i*5.8)+5.7])(window0) for i in range(7)])
plane01=DIFFERENCE([MKPOL([pts1,cells0,pols0]), windows00])
plane1=DIFFERENCE([plane01, windows01])

n=STRUCT([T([1,2,3])([0,0,53])(PROD([plane, Q(-0.4)])),T([1,2,3])([0,0,46.8])(PROD([plane1, Q(-0.3)]))])
s=STRUCT([T([1,2,3])([0,0,0])(PROD([plane, Q(0.4)])),T([1,2,3])([0,0,6.2])(PROD([plane1, Q(0.3)]))])

w=STRUCT([T([1,2,3])([0,0,53])(PROD([plane, Q(-0.4)])),T([1,2,3])([0,0,46.8])(PROD([plane1, Q(-0.3)]))])
e=STRUCT([T([1,2,3])([0,0,0])(PROD([plane, Q(0.4)])),T([1,2,3])([0,0,6.2])(PROD([plane1, Q(0.3)]))])
nord=MAP([S3,S1,S2])(n)

sud=MAP([S3,S1,S2])(s)

est=MAP([S1,S3,S2])(w)
ovest=MAP([S1,S3,S2])(e)
north=nord
south=sud
east=est
west=ovest

prospects=STRUCT([north,south,east,west])

VIEW(STRUCT([plants,prospects]))