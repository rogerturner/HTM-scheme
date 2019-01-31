#!chezscheme

;; === HTM-scheme Concept Copyright 2019 Roger Turner. ===
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
  ;;   Boolean, Number, Integer, Fixnum, (Listof X), (Vectorof X) ... = Scheme types
  ;;   X, Y, Z      type parameters (arbitrary types in function type specification etc)
  ;;   X Y -> Z     function with argument types X, Y and result type Z
  ;;   {X}          abbreviation for (Listof X)
  ;;   Nat          natural number (including zero) (Scheme Fixnum or exact Integer)
  ;;   Perm         Fixnum 0-999 interpreted as a permanence value 0.0-0.999
  ;;   Synapse      Fixnum combining index of input/presynaptic cell and Perm value
  ;;   Synapses     vector of Synapse; sorted to facilitate finding pre cell indices
  ;;   Segment      record with Synapses plus management fields
  ;;   Connections  vector of Segment(s) indexed by cell number in Layer
  ;;   Layer        record with algorithm parameters and Connections of layer cells 
  ;;   Macrocolumn  structure of Layers with interconnections (cortical column)
  ;;   Patch        multiple Macrocolumns with lateral connections
  ;;
  ;; Perm and Synapse types are intended to minimize space required for synapses.
  ;;
  ;; Algorithms have been translated from numenta /htmresearch/... and /nupic/...
  ;; using corresponding functions, names, and organization (not idiomatic Scheme).
  ;; To review HTM-scheme algorithms view Scheme and corresponding Python side-by-side:
  ;; see comments in Numenta code for descriptions of functions and parameters.
  ;; Libraries generally export plain names: it is expected that using modules will
  ;; prefix or rename when importing. Abbreviations are used for record type names to
  ;; avoid excessively long accessor forms.
  ;; Indentation facilitates using a "Fold All" view (in eg Atom) for an overview.

(library (HTM-scheme HTM-scheme algorithms htm_concept)
                                                                                            ;
(export
  max-perm
  min-perm
  clip-max
  clip-min
  perm
  make-synapses
  synapses-ref
  synapses-set!
  synapses-length
  make-syn
  syn-prex
  syn-perm
  build-synapses
  synapses-search
  make-seg
  seg-cellx
  seg-flatx
  seg-synapses
  seg-synapses-set!
  )
                                                                                            ;
(import 
  (except (chezscheme) add1 make-list random reset)
  (HTM-scheme HTM-scheme algorithms htm_prelude))

;; --- Permanence values ---
                                                                                            ;
(define max-perm (- fx3 1))
(define min-perm 0)
                                                                                            ;
(define (clip-max perm)                  ;; Perm -> Perm
  (fxmin max-perm perm))
                                                                                            ;
(define (clip-min perm)                  ;; Perm -> Perm
  (fxmax min-perm perm))
                                                                                            ;
(define (perm x)                         ;; Number[0.0-.999] -> Perm
  (clip-max (clip-min (fx3<- x))))

;; --- Synapses ---
                                                                                            ;
(define-syntax make-synapses             ;; Nat [Fixnum] -> Synapses
  (identifier-syntax make-fxvector))     ;; (change fxvector to vector for R6RS)
                                                                                            ;
(define-syntax synapses-ref              ;; Synapses Nat -> Synapse
  (identifier-syntax fxvector-ref))
                                                                                            ;
(define-syntax synapses-set!             ;; Synapses Nat Synapse ->
  (identifier-syntax fxvector-set!))
                                                                                            ;
(define-syntax synapses-length           ;; Synapses -> Nat
  (identifier-syntax fxvector-length))
                                                                                            ;
(define prex-shift 1024)
                                                                                            ;
(define prex-offset                      ;; -> Fixnum
  ;; offset synapse values to maximize available pre-cell index
  (add1 (fxdiv (least-fixnum) prex-shift)))
                                                                                            ;
(define (make-syn prex perm)             ;; PreX Perm -> Synapse
  (fx+ (fx* (fx+ prex prex-offset) prex-shift) perm))
                                                                                            ;
(define (syn-prex synapse)               ;; Synapse -> PreX
  (fx- (fxdiv synapse prex-shift) prex-offset))
                                                                                            ;
(define (syn-perm synapse)               ;; Synapse -> Perm
  (fxmod synapse prex-shift))
                                                                                            ;
(define (build-synapses n proc)          ;; Nat (Nat -> Synapse) -> Synapses
  (let ((synapses (make-synapses n)))
    (do ((i 0 (add1 i))) ((fx=? i n) synapses)
      (synapses-set! synapses i (proc i)))))
                                                                                            ;
(define (synapses-search                 ;; Synapses Synapse Synapse -> Synapse|#f
          syns syn-low syn-high)
  (let search ((left 0) (right (fx- (synapses-length syns) 1)))
    (if (fx>? left right) #f
      (let* ( (mid (fxdiv (fx+ left right) 2))
              (synapse (synapses-ref syns mid))) 
        (cond 
          [ (fx<? synapse  syn-low) (search (add1 mid) right) ]
          [ (fx<? syn-high synapse) (search left (fx- mid 1)) ]
          [else synapse])))))

;; --- Segments ---
                                                                                            ;
(define-record-type seg                  ;; Segment
  (fields
    cellx                                ;; Nat: index of cell that this is a segment of
    flatx                                ;; Nat: index of this in the register of segments
    (mutable synapses)))                 ;; (vectorof Synapse)

)