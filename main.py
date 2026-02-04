import Optics
import matplotlib.pyplot as plt 

k = (1,-1)

rayon = Optics.Laser(633, k, stokes_vector=[1,1,0,0])

table = Optics.TableOptique((20,20))

mirror1 = Optics.Mirror((1,0.01), 1)
mirror2 = Optics.Mirror((1,1), 1)
mirror3 = Optics.Mirror((1,1), 1)
mirror4 = Optics.Mirror((1,10), 1)
mirror5 = Optics.Mirror((0.5,0.5), 1)
mirror6 = Optics.Mirror((1,0.01), 1)

BS1 = Optics.BeamSplitter((0,1),ratio=0.3)
BS2 = Optics.BeamSplitter((0,1),ratio=0.3)
BS3 = Optics.BeamSplitter((-1,1),ratio=0.3)

P2 = Optics.Polarizer((1,0),90) 
P4 = Optics.Polarizer((1,0),-45) 

P3 = Optics.Mirror((0.25,1),1) 

# table.add(mirror1, (7,0))
# table.add(mirror2, (1.8,-6.5))
# table.add(mirror4, (8.55,-2.2))
table.add(mirror5, (4,0))
# table.add(mirror6, (1,6))

table.add(BS1,(0,-5))
table.add(mirror1,(2,-2))

# table.add(BS2,(5,-5))
# table.add(BS3,(6,-6))

# table.add(BS1,(4,-4))

# table.add(P2,(1.5,-5))
# table.add(P3,(1.7,-8))
# table.add(P4,(5.79, -4.79))

table.add(P2, (5 ,0))
# table.add(P4, (6, 0))
# table.remove(P2)
# table.add(P4, (0, -2))
# table.add(mirror6, (-2.5, -2.5))
# table.add(P2, (-1.41, 4.68))
# table.add(P2, (-2.5, -2.5))

table.draw()
table.draw_laser(rayon, (0,0))
table.report(rayon, (0,0), "test.txt")
plt.show()