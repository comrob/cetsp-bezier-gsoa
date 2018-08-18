/*
 * Filename: ring.cc
 * Author:   Jan Faigl, Petr Vana
 * Date:     2017/06/11 15:34
 */

#include <complex>

#include <crl/logging.h>

#include "neuron.h"
#include "ring.h"

using crl::logger;

using namespace gsoa;

/// - constructor --------------------------------------------------------------
SRing::SRing(void) : neurons(0), nbrNeurons(0) {}

/// - destructor ---------------------------------------------------------------
SRing::~SRing() { deallocate_neurons(); }

/// - public method ------------------------------------------------------------
void SRing::allocate_neurons(int m)
{
   deallocate_neurons();
   neurons = new SNeuron();
   neurons->next = neurons;
   neurons->prev = neurons;
   nbrNeurons = 1;
   for (int i = 1; i < m; i++) {
      insertNeuron(neurons, new SNeuron());
   }
}

/// - public method ------------------------------------------------------------
void SRing::deallocate_neurons(void)
{
   while (nbrNeurons > 0) {
      removeNeuron(neurons->next);
   }
}

/// - public method ------------------------------------------------------------
void SRing::removeNeuron(SNeuron *neuron)
{
   neuron->prev->next = neuron->next;
   neuron->next->prev = neuron->prev;
   if (neurons == neuron) {
      neurons = neuron->next;
      if (nbrNeurons == 1) {
         neurons = 0;
      }
   }
   delete neuron;
   nbrNeurons--;
}

/// - public method ------------------------------------------------------------
SNeuron *SRing::insertNeuronBefore(SNeuron *cur, SNeuron *neuron)
{
   return insertNeuron(cur->next, neuron);
}

/// - public method ------------------------------------------------------------
SNeuron *SRing::insertNeuron(SNeuron *cur, SNeuron *neuron)
{
   // add neuron after the cur neuron
   neuron->prev = cur;
   neuron->next = cur->next;
   cur->next = neuron;
   neuron->next->prev = neuron;
   nbrNeurons++;
   return neuron;
}

/// - public method ------------------------------------------------------------
void SRing::init(const CoordsVector &pts)
{
   SNeuron *n = neurons;
   for (int i = 0; i < pts.size(); ++i) {
      n->setPoint(pts[i]);
      n = n->next;
   }
}

/// - public method ------------------------------------------------------------
void SRing::reset(void)
{
   SNeuron *n = neurons;
   for (int i = 0; i < nbrNeurons; i++) {
      n = n->next;
   }
}

/* end of ring.cc */
