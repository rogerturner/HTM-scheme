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
  ;;   Perm         Fixnum 0-255 interpreted as a permanence value 0.0-1.0
  ;;   Synapse      Fixnum combining index of input/presynaptic cell and Perm value
  ;;   Synapses     Bytevector of Synapse; 32-bit elements; sorted
  ;;   Segment      Record with Synapses plus management fields
  ;;   Connections  Vector of {Segment} indexed by cell number in Layer
  ;;   Layer        Record with algorithm parameters and Connections of layer cells 
  ;;   Macrocolumn  structure of Layers with interconnections (cortical column)
  ;;   Patch        multiple Macrocolumns with lateral connections
  ;;
  ;; Perm and Synapse types are intended to minimize space required for synapses.
  ;; Synapses are sorted for finding pre cell indices; elements of "pre-index" vectors
  ;; (indexed by pre-cell/input) are lists of segments containing the synapse.
  ;;
  ;; Algorithms have been translated from numenta /htmresearch/... and /nupic/...
  ;; using corresponding functions, names, and organization (not idiomatic Scheme).
  ;; To review HTM-scheme algorithms view Scheme and corresponding Python side-by-side:
  ;; see comments in Numenta code for descriptions of functions and parameters.
  ;; Libraries generally export plain names: it is expected that using modules will
  ;; prefix or rename when importing. Abbreviations are used for record type names to
  ;; avoid long accessor forms.
  ;;
  ;; Function definitions are annotated with their type and a brief description;
  ;; indentation facilitates using a "Fold All" view (in eg Atom) for an overview.
  ;;
  ;; "Lisp's parentheses are the bumps on the top of Lego" [Paul Graham]
  ;; "The only honest function names are car and cdr"  [citation needed]

(library (HTM-scheme HTM-scheme algorithms htm_concept)
                                                                                            ;
(export
  max-perm
  min-perm
  clip-max
  clip-min
  make-syn
  syn-prex
  syn-perm
  perm
  make-seg
  seg-cellx
  seg-overlap
  seg-overlap-set!
  seg-overlap-clear!
  seg-synapses
  seg-synapses-set!
  make-synapses
  synapses-ref
  synapses-set!
  synapses-length
  seg-synapses->list
  seg-synapses<-list
  synapses-search
  expect
  )
                                                                                            ;
(import 
  (chezscheme)
  (HTM-scheme HTM-scheme algorithms htm_prelude))

  ;; Convenience abbreviations
  (define-syntax fxasl  (identifier-syntax fxarithmetic-shift-left))
  (define-syntax fxasr  (identifier-syntax fxarithmetic-shift-right))
  (define-syntax bvu32  (identifier-syntax bytevector-u32-native-ref))
  (define-syntax bvu32! (identifier-syntax bytevector-u32-native-set!))

;; --- Permanence and Synapse values ---
                                                                                            ;
(define prex-shift 8)
                                                                                            ;
(define min-perm   0)
                                                                                            ;
(define max-perm
  (- (expt 2 prex-shift) 1))
                                                                                            ;
(define (make-syn prex perm)             ;; PreX Perm -> Synapse
  (fx+ (fxasl prex prex-shift) perm))
                                                                                            ;
(define (syn-prex synapse)               ;; Synapse -> PreX
  (fxasr synapse prex-shift))
                                                                                            ;
(define (syn-perm synapse)               ;; Synapse -> Perm
  (fxand synapse max-perm))
                                                                                            ;
(define (clip-max perm)                  ;; Perm -> Perm
  (fxmin max-perm perm))
                                                                                            ;
(define (clip-min perm)                  ;; Perm -> Perm
  (fxmax min-perm perm))
                                                                                            ;
(define (perm x)                         ;; Number[0.0-1.0] -> Perm
  (if (zero? x)  0
    (clip-max (fxmax 1 (int<- (* x max-perm))))))

;; --- Segments ---
                                                                                            ;
(define-record-type seg
  (fields
    cellx                         ;; Nat: index of cell that this is a segment of
    (mutable overlap)             ;; synapses bv length-4; overlap counts
    (mutable synapses))           ;; Bytevector: 32-bit synapses
  (protocol
    (lambda (new)
      (lambda (cellx)
        (new cellx (fxasl -4 32) (make-synapses))))))
                                                                                            ;
(define (seg-overlap-clear! s)
  ;; clear overlap counts only, preserve synapses length
  (seg-overlap-set! s (fxand #xFFFF00000000 (seg-overlap s))))
                                                                                            ;
(define (make-synapses)
  (make-bytevector 0))
                                                                                            ;
(define (synapses-ref bv n)
  (bvu32 bv (fx* n 4)))
                                                                                            ;
(define (synapses-set! bv n u32)
  (bvu32! bv (fx* n 4) u32))
                                                                                            ;
(define (synapses-length bv)
  (fxdiv (bytevector-length bv) 4))
                                                                                            ;
(define native (native-endianness))
                                                                                            ;
(define (seg-synapses->list seg)         ;; Segment -> {Synapse}
  (bytevector->uint-list (seg-synapses seg) native 4))
                                                                                            ;
(define (seg-synapses<-list seg ss)      ;; Segment {Synapse} ->
  ;; replace seg's synapses with ss (overlap counts are not preserved)
  (let ([synapses (uint-list->bytevector ss native 4)])
    (seg-overlap-set! seg (fxasl (fx- (bytevector-length synapses) 4) 32))
    (seg-synapses-set! seg synapses)))
                                                                                            ;
(define-syntax synapses-search           ;; Segment PreX -> Synapse | #f
  ;; binary search for synapse (forced inline)
  (lambda (x)
    (syntax-case x ()
      [ (_ seg prex)
        #'(let* ( [limit  (fxasr (seg-overlap seg) 32)]
                  [syns   (seg-synapses seg)]
                  [target (fxasl prex prex-shift)])
            (let search ([left 0] [right limit])
              (if (fx<=? left right)
                (let* ( [mid      (fxand #xFFFFFC (fxasr (fx+ left right) 1))]
                        [synapse  (bvu32 syns mid)])
                  (cond 
                    [ (fx>? target synapse) (search (fx+ mid 4) right)]
                    [ (fx<? target (fxand #xFFFFFF00 synapse)) (search left (fx- mid 4)) ]
                    [else synapse]))
                #f))) ])))

(define-syntax expect                    ;; ((X ... -> Y) X ...) Y -> Error|
  ;; check that function application(s) to arguments match expected values
  (lambda (x)                            
    (syntax-case x ()                    ;; [expect ([fn args] expected ) ... ]
      [ (_ (expr expected) ...)          ;; expr matches [fn args]
        #'(begin 
            (let ((result expr))         ;; eval expr just once, no output if check passes
              ; (when (equal? result expected) (display "."))
              (unless (equal? result expected)
                (for-each display 
                  `("**" expr #\newline 
                    "  expected: " ,expected #\newline 
                    "  returned: " ,result  #\newline))
                (exit))) ...)])))

)