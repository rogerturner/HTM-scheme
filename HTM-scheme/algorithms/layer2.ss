#!chezscheme

;; === HTM-scheme Layer2 algorithm  (C) 2021 Roger Turner. ===
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

  ;; See htm_concept.ss for type and data structure description and code conventions.
  
(library (HTM-scheme HTM-scheme algorithms layer2)
                                                                                            ;
(export
make-l2
compute
reset
get-active-cells
get-active-cols
print-statistics)
                                                                                            ;
(import
  (except (chezscheme) reset)
          (parameters)
          (HTM-scheme HTM-scheme algorithms htm_prelude)
          (HTM-scheme HTM-scheme algorithms htm_concept)
  (prefix (HTM-scheme HTM-scheme algorithms column_pooler) cp:))
                                                                                            ;
  (implicit-exports #f)

;; layer coding in Source (cells that project to layer2):
;;    1  p23
;;    3  ss4l23
;;    4  p4

;; === Layer record ===
                                                                                            ;
(define-record-type l2 (fields           ;; L2
    p23)                                 ;; CP
(protocol #;(make-l2 ncc ccx nmc l2-c/mc ss4l23-c/mc p4-c/mc cp-overrides)
  ;; configure cp instance with connect? proc
  (lambda (new)
    (lambda (
        ncc                              ;; number of cortical columns
        ccx                              ;; index of this cortical column
        nmc                              ;; number of minicolumns/cc
        l2-c/mc                          ;; cells/minicolumn
        cp-overrides connect?)
      (new
        (cp:make-cp (append `(
                        [seed                 . ,(random 4294967295)]
                        [cell-count           . ,(fx* l2-c/mc nmc)]
                        [cortical-column      . ,ccx]
                        [lateral-input-widths . ,(make-list ncc 0)]
                        [connect?             . ,(lambda (source post-mcx)
                                                   (connect? source post-mcx ccx))])
                      cp-overrides)))))))

;; === Layer algorithm ===
                                                                                            ;
(define  compute                         ;; ...
  ;; run one time step of the column pooler algorithm
  (case-lambda
#;> [(l2 feedforward-input learn)        ;; L2 {CellX} Boolean ->
      (compute l2 feedforward-input '() '() '() learn '() '()) ]
                                                                                            ;
#;> [(l2 feedforward-input distal-input  ;; L2 {CellX} {{CellX}} {CellX} Boolean ->
        lateral-inputs feedforward-growth-candidates learn)
      (compute l2 feedforward-input distal-input lateral-inputs feedforward-growth-candidates learn '() '()) ]
                                                                                            ;
#;> [(l2 feedforward-input distal-input  ;; L2 {CellX} {CellX} {{CellX}} {CellX} Boolean {CellX} {ColX} ->
        lateral-inputs feedforward-growth-candidates learn predicted-input inhibited-cols)
      (cp:compute (l2-p23 l2)
        feedforward-input distal-input lateral-inputs feedforward-growth-candidates
        learn predicted-input inhibited-cols) 
      ]))
                                                                                            ;
(define (reset l2)                       ;; L2 ->
  (cp:reset (l2-p23 l2)))

;; === Accessors ===
                                                                                            ;
(define (get-active-cells l2)
  (cp:get-active-cells (l2-p23 l2)))
                                                                                            ;
(define (get-active-cols l2)
  (cp:get-active-cols (l2-p23 l2)))
                                                                                            ;
(define (print-statistics sum-fl)
  (define (sum-l2 f) (sum-fl f l2-p23))
  (when (positive? (sum-l2 cp:cp-n-sdrs))
    (for-each display `(
        ,(sum-l2 cp:cp-n-sdrs) " sdrs created\nsyns/segments   "
        ,(sum-l2 cp:number-of-proximal-synapses) "/"
        ,(sum-l2 cp:number-of-proximal-segments) " p23 proximal\n"))
    (when (positive? (sum-l2 cp:number-of-distal-synapses))
      (for-each display `(
          "                " ,(sum-l2 cp:number-of-distal-synapses) "/"
                             ,(sum-l2 cp:number-of-distal-segments) " p23 distal\n")))
    (when (positive? (sum-l2 cp:number-of-lateral-synapses))
      (for-each display `(
          "                " ,(sum-l2 cp:number-of-lateral-synapses) "/"
                             ,(sum-l2 cp:number-of-lateral-segments) " p23 lateral\n")))))

)