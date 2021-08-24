#!/usr/bin/python
import argparse
import numpy as np
import random
import itertools
import math
import scipy.signal
import time
import matplotlib.pyplot as plt

from grid_cells import GridCells
from nupic.encoders.coordinate import CoordinateEncoder


class Environment(object):
  """
  Environment is a 2D square, in first quadrant with corner at origin.
  """
  def __init__(self, size):
    self.size     = size
    self.position = (size/2, size/2)
    self.speed    = 2.0 ** .5
    self.angle    = 0
    self.course   = []

  def in_bounds(self, position):
    x, y = position
    x_in = x >= 0 and x < self.size
    y_in = y >= 0 and y < self.size
    return x_in and y_in

  def move(self):
    max_rotation = 2 * math.pi / 20
    self.angle += random.uniform(-max_rotation, max_rotation)
    vx = self.speed * math.cos(self.angle)
    vy = self.speed * math.sin(self.angle)
    x, y = self.position
    new_position = (x + vx, y + vy)
    
    if self.in_bounds(new_position):
      self.position = new_position
      self.course.append(self.position)
    else:
      # On failure, recurse and try again.
      assert(self.in_bounds(self.position))
      self.angle = random.uniform(0, 2 * math.pi)
      self.move()

  def plot_course(self, show=True):
    plt.figure("Path")
    plt.ylim([0, self.size])
    plt.xlim([0, self.size])
    x, y = zip(*self.course)
    plt.plot(x, y, 'k-')
    if show:
      plt.show()


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--train_time', type=int, default = 1000 * 1000,)
  parser.add_argument('--debug', action='store_true')
  args = parser.parse_args()

  # Setup
  env = Environment(size = 200)
  enc = CoordinateEncoder(w=75, n=2500)
  assert(enc.w > 50)
  enc_radius = 5
  gcm = GridCells(
    b1               = 0.05,
    b2               = 0.05 / 3,
    inputDimensions  = (enc.n,),
    columnDimensions = (100,),
    potentialPct     = 0.95,
    numActiveColumnsPerInhArea = .2 * 100,
    synPermInactiveDec = 0.0008,
    synPermActiveInc   = 0.005,
    synPermConnected   = 0.25,
    stimulusThreshold  = 0,
    boostStrength      = 0.,
    globalInhibition = True,
    potentialRadius  = enc.n,
    wrapAround       = True,)


  def compute(learn=True):
    nearest_position = np.array(np.rint(env.position), dtype=np.int)
    enc_sdr = np.zeros(enc.n)
    enc.encodeIntoArray((nearest_position, enc_radius), enc_sdr)

    gc_act = np.zeros(gcm.getNumColumns())
    gcm.compute(enc_sdr, learn, gc_act)
    return enc_sdr, gc_act


  print("Training for %d cycles ..."%args.train_time)
  start_time = time.time()
  gcm.reset()
  for step in range(args.train_time):
    if step % 1000 == 0:
      print("Cycle %d"%step)
    env.move()
    compute()
  train_time = time.time()
  print("Elapsed time (training): %d seconds."%int(round(train_time - start_time)))

  print("Testing ...")

  # Show how the agent traversed the environment.
  env.plot_course(show=False)

  # Measure Receptive Fields.
  enc_num_samples = 12
  gc_num_samples  = 20
  enc_samples = random.sample(xrange(enc.n), enc_num_samples)
  gc_samples  = random.sample(xrange(gcm.getNumColumns()), gc_num_samples)
  enc_rfs = [np.zeros((env.size, env.size)) for idx in enc_samples]
  gc_rfs  = [np.zeros((env.size, env.size)) for idx in gc_samples]
  for position in itertools.product(xrange(env.size), xrange(env.size)):
    env.position = position
    gcm.reset()
    enc_sdr, gc_sdr = compute(learn=False)
    for rf_idx, enc_idx in enumerate(enc_samples):
      enc_rfs[rf_idx][position] = enc_sdr[enc_idx]
    for rf_idx, gc_idx in enumerate(gc_samples):
      gc_rfs[rf_idx][position] = gc_sdr[gc_idx]
  
  # Show the Input/Encoder Receptive Fields.
  if enc_num_samples > 0:
    plt.figure("Input Receptive Fields")
    nrows = int(enc_num_samples ** .5)
    ncols = math.ceil((enc_num_samples+.0) / nrows)
    for subplot_idx, rf in enumerate(enc_rfs):
      plt.subplot(nrows, ncols, subplot_idx + 1)
      plt.imshow(rf)

  # Show the Grid Cells Receptive Fields.
  if gc_num_samples > 0:
    plt.figure("Grid Cell Receptive Fields")
    nrows = int(gc_num_samples ** .5)
    ncols = math.ceil((gc_num_samples+.0) / nrows)
    for subplot_idx, rf in enumerate(gc_rfs):
      plt.subplot(nrows, ncols, subplot_idx + 1)
      plt.imshow(rf)

    # Show the autocorrelations of the grid cell receptive fields.
    plt.figure("Grid Cell RF Autocorrelations")
    for subplot_idx, rf in enumerate(gc_rfs):
      plt.subplot(nrows, ncols, subplot_idx + 1)
      xcor = scipy.signal.correlate2d(rf, rf)
      plt.imshow(xcor)

  test_time = time.time()
  print("Elapsed time (testing): %d seconds."%int(round(test_time - train_time)))
  plt.show()
