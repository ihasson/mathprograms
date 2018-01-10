
(define (fixed-point-iter fnct p0 TOL maxiter)  
  (define (fpi p i)
    (let ([pn (fnct p)])
      (cond 
        [(= i maxiter) (displayln 
                         (list "Max iterations " maxiter  " reached. "
                               "Last value " pn)) ]
        [(< (abs (- pn  p)) TOL) (begin (displayln pn) pn)]
        [else (fpi pn (+ i 1))])))
  (fpi p0 0))

(define (g5 x) (- x (/ [+ (* x x x) (* 4 x x) -10] [+ (* 3 x x) (* 8 x) ] ))) 
(define (g1 x) (+ x (* -1 x x x) (* -4 x x) 10))
(define (g2 x) (sqrt (+ (/ 10 x) (* -4 x))))
(define (g3 x) (* 0.5 (sqrt (- 10 (expt x 3)))))
(define (g4 x) (sqrt (/ 10 (+ x 4))))

