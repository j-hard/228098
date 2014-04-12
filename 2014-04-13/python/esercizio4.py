from pyplasm import *
from math import sqrt  
import sys
sys.path.insert(0, '/home/j-hard/lar-cc/lib/py/')
from simplexn import *
from larcc import *
from lar2psm import * 
from largrid import *
from scipy import *

pts0=[[-20,-20],[73,-20],[73,73],[-20,73]]
ptsbasestrada=[[-200,-300],[400,-300],[400,300],[-200,300]]
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
pts3=[[5.8,0],[0.8,0],[0.8,6.6],[5.8,6.6]]
pts4=[[5.8,0],[0.8,0],[0.8,5.6],[5.8,5.6]]

def disk2D(r):
	u,v=r
	return [v*COS(u),v*SIN(u)]


d2d=PROD([INTERVALS(2*PI)(8),INTERVALS(2.5)(3)])
disk=MAP(disk2D)(d2d)

window0=STRUCT([MKPOL([pts3,cells0,pols0]),T([1,2])([3.3,6.6])(disk)])
windows0=STRUCT([T([1])([(i*5.8)])(window0) for i in range(9)])
plane0=DIFFERENCE([MKPOL([pts0,cells0,pols0]), windows0])

windowl1=STRUCT([MKPOL([pts4,cells0,pols0]),T([1,2])([3.3,5.6])(disk)])
windowsl1=STRUCT([T([1,2])([(i*5.8),10])(windowl1) for i in range(9)])
windows1=STRUCT([T([2])([(i*9)])(windowsl1) for i in range(5)])
plane=DIFFERENCE([plane0, windows1])

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

pts0=[[-20,-20],[70,-20],[70,70],[-20,70]]
pts1=[[0,0],[53,0],[53,53],[0,53]]
pts2=[[6.2,6.2],[46.8,6.2],[46.8,46.8],[6.2,46.8]]
cells=[[1,2,3,4]]
pols=None


ptsscala=[[0,0],[0.4,0],[0.4, -0.2],[0,-0.2]]
ptscorrimano=[[0,0.2],[12,0.2],[12, -5.8],[0,-5.8]]
muroperimetrale=[[0,0.2],[93,0.2],[93, -5.8],[0,-5.8]]

gradino3D=PROD([MKPOL([ptsscala,cells,pols]), Q(93)])
corrimano3D=PROD([MKPOL([ptscorrimano,cells,pols]), Q(1)])
muroperimetrale3D=PROD([MKPOL([muroperimetrale,cells,pols]), Q(0.4)])

scala3D0=STRUCT([T([1,2,3])([0.4*(x-1),-0.2*(x-1),0])(gradino3D) for x in range(30)])
scala3D=STRUCT([scala3D0,corrimano3D,T([1,2,3])([0,0,92])(corrimano3D),T([1,2,3])([-93,0,0])(muroperimetrale3D)])
scalaEst=T([1,2,3])([-20,73,-0.2])(MAP([S3,S1,S2])(scala3D))

scalaNord=T([1,2,3])([0,53,0])(R([1,2])(-PI/2)(scalaEst))
scalaSud=T([1,2,3])([53,0,0])(R([1,2])(PI/2)(scalaEst))

scale=STRUCT([scalaNord,scalaSud])
model3D=STRUCT([plants,prospects,scale])


palazzo1=CUBOID([230,30,25])
palazzo2=CUBOID([100,30,15])
palazzo3=CUBOID([30,100,25])
palazzo4=CUBOID([50,20,15])
palazzo5=CUBOID([30,30,40])
palazzo6=CUBOID([40,110,35])
pal7a=CUBOID([45,15,15])
pal7b=CUBOID([15,15,15])
palazzo7=STRUCT([pal7a,T([1,2,3])([15,-15])(pal7b)])
pal8a=CUBOID([20,15,15])
pal8b=CUBOID([15,sqrt(450),15])
palazzo8=STRUCT([pal8a,R([1,2])(PI/4)(pal8b),T([1,2,3])([0,15,0])(R([1,2])(PI/2)(pal8a))])

basestrada=MKPOL([ptsbasestrada,cells,pols])

davanti=STRUCT([T([1,2,3])([85,-100,0])(palazzo1),T([1,2,3])([135,-50,0])(palazzo2),T([1,2,3])([135,73,0])(palazzo2)])
lato=STRUCT([T([1,2,3])([35,-160,0])(palazzo4),T([1,2,3])([35,-250,0])(palazzo4),T([1,2,3])([65,-210,0])(palazzo5),T([1,2,3])([-62,-170,0])(palazzo3),T([1,2,3])([225,-250,0])(palazzo6)])
gruppo=STRUCT([T([1,2,3])([115,-155,0])(palazzo7),T([1,2,3])([115,-235,0])(palazzo7),T([1,2,3])([195,-185,0])(R([1,2])(PI/2)(palazzo7)),T([1,2,3])([210,-235,0])(R([1,2])(PI/2)(palazzo8))])
edifici=STRUCT([T([1,2,3])([0,0,6])(model3D),davanti,lato,gruppo,COLOR(BLACK)(basestrada)])

ptspianta0=[[-20,0],[20,0],[20,-100],[-20,-100]]
ptspianta1=[[0,0],[10,-10],[0,-20],[-10,-10]]
ptspianta1b=[[0,-2],[8,-10],[0,-18],[-8,-10]]
ptspianta2=[[0,0],[5,-5],[0,-10],[-5,-5]]
ptspianta2b=[[0,-2],[3,-5],[0,-8],[-3,-5]]
basepianta0=MKPOL([ptspianta0,cells,pols])
basepianta1=MKPOL([ptspianta1,cells,pols])
basepianta1b=MKPOL([ptspianta1b,cells,pols])
basepianta2=MKPOL([ptspianta2,cells,pols])
basepianta2b=MKPOL([ptspianta2b,cells,pols])
pianta0=PROD([basepianta0, Q(0.1)])
pianta1=DIFFERENCE([PROD([basepianta1, Q(1)]),PROD([basepianta1b, Q(1)])])
pianta2=DIFFERENCE([PROD([basepianta2, Q(1)]),PROD([basepianta2b, Q(1)])])
aiuola=STRUCT([pianta0,T([1,2])([0,-10])(pianta1),T([1,2])([0,-40])(pianta1),T([1,2])([0,-70])(pianta1),pianta2,T([1,2])([0,-30])(pianta2),T([1,2])([0,-60])(pianta2),T([1,2])([0,-90])(pianta2)])

ptsbasealbero=[[-1,1],[1,1],[1,-1],[-1,-1]]
chiomalbero=T([1,2,3])([-2.5,-2.5,-2.5])(CUBOID([5,5,5]))
basealbero=MKPOL([ptsbasealbero,cells,pols])
sezionetronco=PROD([INTERVALS(2*PI)(8),INTERVALS(0.5)(3)])
tronco=PROD([MAP(disk2D)(sezionetronco), Q(10)])
chioma=STRUCT([chiomalbero,R([1,2])(PI/4)(chiomalbero),R([1,3])(PI/4)(chiomalbero),R([2,3])(PI/4)(chiomalbero)])
albero=STRUCT([COLOR(BROWN)(tronco),COLOR(GREEN)(basealbero),COLOR(GREEN)(T([1,2,3])([0,0,10])(chioma))])
filalberi=STRUCT([T([1,2])([((i)*10),0])(albero) for i in range(10)])
filalberi1=STRUCT([T([1,2])([((i)*10),0])(albero) for i in range(5)])
gruppoalberi=STRUCT([filalberi1,T([1,2,3])([0,10,0])(filalberi1),T([1,2,3])([0,20,0])(filalberi1)])
vegetazione=STRUCT([COLOR(GREEN)(T([1,2])([135,26.5])(R([1,2])(PI/2)(aiuola))),T([1,2,3])([137.5,-10,0])(filalberi),T([1,2,3])([137.5,63,0])(filalberi),T([1,2,3])([-22,-167.5,0])(R([1,2])(PI/2)(filalberi)),T([1,2,3])([25,-247.5,0])(R([1,2])(PI/2)(filalberi)),T([1,2,3])([25,-100,0])(gruppoalberi),T([1,2,3])([115,-205,0])(gruppoalberi)])
quartiere=STRUCT([edifici,vegetazione])

VIEW(quartiere)