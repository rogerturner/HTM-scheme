#!/usr/bin/python

import numpy
from nupic.algorithms.spatial_pooler import SpatialPooler, realDType
from nupic.algorithms.temporal_memory_shim import TemporalMemoryShim as TemporalMemory

class GridCells(SpatialPooler):
  """
  """
  def __init__(self, b1, b2, *args, **kw_args):
    self.b1 = b1
    self.b2 = b2
    SpatialPooler.__init__(self, *args, **kw_args)
    self.reset()

  def reset(self):
    self.R_act           = numpy.zeros(self.getNumColumns(), dtype=realDType)
    self.R_inact         = numpy.zeros(self.getNumColumns(), dtype=realDType)
    self.prevPermChanges = {}

  def compute(self, inputVector, learn, activeArray):
    """
    This is the primary public method of the SpatialPooler class. This
    function takes a input vector and outputs the indices of the active columns.
    If 'learn' is set to True, this method also updates the permanences of the
    columns.

    :param inputVector: A numpy array of 0's and 1's that comprises the input
        to the spatial pooler. The array will be treated as a one dimensional
        array, therefore the dimensions of the array do not have to match the
        exact dimensions specified in the class constructor. In fact, even a
        list would suffice. The number of input bits in the vector must,
        however, match the number of bits specified by the call to the
        constructor. Therefore there must be a '0' or '1' in the array for
        every input bit.
    :param learn: A boolean value indicating whether learning should be
        performed. Learning entails updating the  permanence values of the
        synapses, and hence modifying the 'state' of the model. Setting
        learning to 'off' freezes the SP and has many uses. For example, you
        might want to feed in various inputs and examine the resulting SDR's.
    :param activeArray: An array whose size is equal to the number of columns.
        Before the function returns this array will be populated with 1's at
        the indices of the active columns, and 0's everywhere else.
    """
    if not isinstance(inputVector, numpy.ndarray):
      raise TypeError("Input vector must be a numpy array, not %s" %
                      str(type(inputVector)))

    if inputVector.size != self._numInputs:
      raise ValueError(
          "Input vector dimensions don't match. Expecting %s but got %s" % (
              inputVector.size, self._numInputs))

    self._updateBookeepingVars(learn)
    inputVector = numpy.array(inputVector, dtype=realDType)
    inputVector.reshape(-1)
    self._overlaps = self._calculateOverlap(inputVector)

    # Divide by the total synaptic input strength.
    sum_synapses = numpy.empty(self.getNumColumns())
    self.getConnectedCounts(sum_synapses)
    self._overlaps = self._overlaps / sum_synapses

    # Apply fatigue and update the fatigue variables
    self._overlaps = self._fatigue(self._overlaps)

    # Apply boosting when learning is on
    if learn:
      self._boostedOverlaps = self._boostFactors * self._overlaps
    else:
      self._boostedOverlaps = self._overlaps

    # Apply inhibition to determine the winning columns
    activeColumns = self._inhibitColumns(self._boostedOverlaps)
    activeColumns = numpy.array(activeColumns, dtype=numpy.int)

    if learn:
      self._adaptSynapses(inputVector, activeColumns)
      self._updateDutyCycles(self._overlaps, activeColumns)
      self._bumpUpWeakColumns()
      self._updateBoostFactors()
      if self._isUpdateRound():
        self._updateInhibitionRadius()
        self._updateMinDutyCycles()

    activeArray.fill(0)
    activeArray[activeColumns] = 1

  def _fatigue(self, overlaps):
    self.R_act   += self.b1 * (overlaps - self.R_inact - self.R_act)
    self.R_inact += self.b2 * (overlaps - self.R_inact)
    return self.R_act

  def _adaptSynapses(self, inputVector, activeColumns):
    """
    The primary method in charge of learning. Adapts the permanence values of
    the synapses based on the input vector, and the chosen columns after
    inhibition round. Permanence values are increased for synapses connected to
    input bits that are turned on, and decreased for synapses connected to
    inputs bits that are turned off.

    Parameters:
    ----------------------------
    :param inputVector:
                    A numpy array of 0's and 1's that comprises the input to
                    the spatial pooler. There exists an entry in the array
                    for every input bit.
    :param activeColumns:
                    An array containing the indices of the columns that
                    survived inhibition.
    """
    inputIndices = numpy.where(inputVector > 0)[0]
    permChanges = numpy.zeros(self._numInputs, dtype=realDType)
    permChanges.fill(-1 * self._synPermInactiveDec)
    permChanges[inputIndices] = self._synPermActiveInc
    allPermChanges = {}
    for columnIndex in activeColumns:
      perm = self._permanences[columnIndex]
      maskPotential = numpy.where(self._potentialPools[columnIndex] > 0)[0]
      columnPermChanges = permChanges[maskPotential]
      allPermChanges[columnIndex] = numpy.array(columnPermChanges, copy=True)
      # Filter out permanence changes which are repeated on consequtive cycles.
      try:
        columnPrevPermChanges = self.prevPermChanges[columnIndex]
      except KeyError:
        columnPrevPermChanges = numpy.zeros(len(maskPotential), dtype=realDType)
      columnPermChanges[columnPermChanges == columnPrevPermChanges] = 0.
      perm[maskPotential] += columnPermChanges
      self._updatePermanencesForColumn(perm, columnIndex, raisePerm=True)
    self.prevPermChanges = allPermChanges
