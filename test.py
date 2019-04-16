import mayavi.mlab as mlab
import numpy as np

Cpos = np.random.random((8, 3, 1700))
Hpos = np.random.random((16, 3, 1700))
allActPos = np.zeros((24, 3, 1700))
qm = np.random.random((24, 3, 1700))

sizes = {'C': 1.8, 'H': 0.7}
pts={}
draw={'pos': True, 'qm': True}

# Plot first step
pts['Cpos'] = mlab.points3d(Cpos[:, 0, 0], Cpos[:, 1, 0], Cpos[:, 2, 0],
                            scale_factor=sizes['C'], color=(0, 0, 0))
pts['Hpos'] = mlab.points3d(Hpos[:, 0, 0], Hpos[:, 1, 0], Hpos[:, 2, 0],
                            scale_factor=sizes['H'], color=(1, 1, 0))
pts['qm'] = mlab.quiver3d(allActPos[:, 0, 0], allActPos[:, 1, 0], allActPos[:, 2, 0],
                          qm[:, 0, 0], qm[:, 1, 0], qm[:, 2, 0])

# Plot animation
@mlab.animate(delay=10)
def anim(draw, pts):
   for step in range(1, Cpos.shape[2]):
      if draw['pos']:
         pts['Cpos'].mlab_source.set(x=Cpos[:, 0, step],
                                    y=Cpos[:, 1, step],
                                    z=Cpos[:, 2, step])
         pts['Hpos'].mlab_source.set(x=Hpos[:, 0, step],
                                    y=Hpos[:, 1, step],
                                    z=Hpos[:, 2, step])
      if draw['qm']:
         pts['qm'].mlab_source.set(x=allActPos[:, 0, step],
                                    y=allActPos[:, 1, step],
                                    z=allActPos[:, 2, step],
                                    u=qm[:, 0, step],
                                    v=qm[:, 1, step],
                                    w=qm[:, 2, step])

      yield

anim(draw, pts)
mlab.show()
