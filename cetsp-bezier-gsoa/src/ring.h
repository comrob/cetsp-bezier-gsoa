/*
 * Filename: ring.h
 * Date:     2017/06/11 15:34
 * Author:   Jan Faigl, Petr Vana
 */

#ifndef __RING_H__
#define __RING_H__

#include "coords.h"

namespace gsoa {

   struct SNeuron;
   struct SRing {
      SNeuron *neurons;
      int nbrNeurons;
      double length;

      SRing(void);
      ~SRing();

      void allocate_neurons(int m);

      void deallocate_neurons(void);

      void removeNeuron(SNeuron *neuron);

      SNeuron *insertNeuronBefore(SNeuron *cur, SNeuron *neuron);

      SNeuron *insertNeuron(SNeuron *cur, SNeuron *neuron);

      void init(const CoordsVector &pts);

      void reset(void);
   };

} // end namespace gsoa

#endif

/* end of ring.h */
