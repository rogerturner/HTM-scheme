#!r6rs

;; === HTM-scheme Concept Copyright 2018 Roger Turner. ===
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; This program is free software: you can redistribute it and/or modify  ;;
  ;; it under the terms of the GNU Affero Public License version 3 as      ;;
  ;; published by the Free Software Foundation.                            ;;
  ;;                                                                       ;;
  ;; This program is distributed in the hope that it will be useful,       ;;
  ;; but WITHOUT ANY WARRANTY; without even the implied warranty of        ;;
  ;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  ;;
  ;; See the GNU Affero Public License for more details.                   ;;
  ;;                                                                       ;;
  ;; You should have received a copy of the GNU Affero Public License      ;;
  ;; along with this program.  If not, see http://www.gnu.org/licenses.    ;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Types and data structure concepts for HTM-scheme algorithms:
  ;;   Perm:        Fixnum 0-9999 interpreted as a permanence value 0.0-0.9999
  ;;   Synapse:     Fixnum containing index of presynaptic cell and permanence value
  ;;   Synapses:    vector of Synapse; sorted to facilitate finding pre cell indices
  ;;   Segment:     record with Synapses plus management fields
  ;;   Connections: vector of Segment(s) indexed by cell number in Layer
  ;;   Layer:       record with algorithm parameters and Connections of layer cells 
  ;;   Macrocolumn: structure of Layers with interconnections
  ;;   Patch:       vector of Macrocolumn with lateral connections
  ;;
  ;; Perm/Synapse/Segment types are intended to minimize space required for synapses.
  ;;
  ;; Algorithms have been translated from numenta /htmresearch/... and /nupic/...
  ;; using corresponding functions, names, and organization (not idiomatic Scheme).
  ;; To review HTM-scheme algorithms view Scheme and corresponding Python side-by-side:
  ;; see comments in Numenta code for descriptions of functions and parameters.
  ;; Libraries generally export plain names: it is expected that using modules will
  ;; prefix or rename when importing. Record type names are short to avoid excessively
  ;; long accessor forms.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(library (HTM-scheme HTM-scheme algorithms htm_concept)
                                                                                            ;
(export
  max-perm
  min-perm
  clip-max
  clip-min
  perm
  syn-least
  make-syn
  syn-prex
  syn-perm
  make-synapses
  synapses-ref
  synapses-set!
  synapses-length
  build-synapses
  make-seg
  seg-cellx
  seg-flatx
  seg-synapses
  seg-synapses-set!
  )
                                                                                            ;
(import 
  (rnrs)
  (HTM-scheme HTM-scheme algorithms htm_prelude))

;; --- Permanence values ---
                                                                                            ;
(define max-perm (- x10k 1))
(define min-perm 0)
                                                                                            ;
(define (clip-max perm)                  ;; Perm -> Perm
  (fxmin max-perm perm))
                                                                                            ;
(define (clip-min perm)                  ;; Perm -> Perm
  (fxmax min-perm perm))
                                                                                            ;
(define (perm x)                         ;; Number[0.0-.9999] -> Perm
  (clip-max (clip-min (int<- (* x x10k)))))

;; --- Synapses ---
                                                                                            ;
(define syn-least (least-fixnum))
                                                                                            ;
(define (make-syn prex perm)             ;; PreX Perm -> Synapse
  (+ (* prex x10k) syn-least perm))
                                                                                            ;
(define (syn-prex synapse)               ;; Synapse -> PreX
  (div (- synapse syn-least) x10k))
                                                                                            ;
(define (syn-perm synapse)               ;; Synapse -> Perm
  (mod (- synapse syn-least) x10k))
                                                                                            ;
(define make-synapses   make-vector)
(define synapses-ref    vector-ref)
(define synapses-set!   vector-set!)
(define synapses-length vector-length)
                                                                                            ;
(define (build-synapses n proc)          ;; Nat (Nat -> Synapse) -> Synapses
  (let ((synapses (make-synapses n)))
    (do ((i 0 (add1 i))) ((fx=? i n) synapses)
      (synapses-set! synapses i (proc i)))))

;; --- Segments ---
                                                                                            ;
(define-record-type seg                  ;; Segment
  (fields
    cellx                                ;; Nat: index of cell that this is a segment of
    flatx                                ;; Nat: index of this in the register of segments
    (mutable synapses)))                 ;; (vectorof Synapse)

)