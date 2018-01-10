;(struct term (deg coef))
(define (term-deg dc)
  (if (pair? dc) (car dc)
    (error "not a term")))
(define (term-coef dc)
  (if (pair? dc) (cadr dc)
    (error "not a term")))


(define (t-sort ls)
  (sort ls (lambda (x y) (> (term-deg x) (term-deg y)))))

(define (t-sum ls acc)
  (cond 
    [(empty? ls) acc]
    [(empty? acc) (t-sum (cdr ls) (cons (car ls) acc))]
    [(= (term-deg (car ls)) (term-deg (car acc)))
     (t-sum (cdr ls) (cons 
                       (list 
                         (term-deg (car acc)) 
                         (+ (term-coef (car acc)) (term-coef (car ls))))
                       (cdr acc)))]
    [else (t-sum (cdr ls) (cons (car ls) acc))]))

(define (reduce ls)
  (let ([s (t-sort ls)])
    (t-sum s '())))

(define (poly-sum p q)
  (reduce (append p q)))

(define (poly-diff p q)
  (reduce (append p (poly-scalar-prod -1 q))))

(define (poly-scalar-prod s p)
  (map (lambda (t) (list (car t) (* s (cadr t)))) p))

(define (poly-prod p q)
  (define (t-poly-prod t pol)
    (map (lambda (pol-t) 
           (list (+ (car t) (car pol-t))
                 (* (cadr t) (cadr pol-t))))
         pol))
  (reduce 
    (foldr append '() 
           (map (lambda (t) (t-poly-prod t q)) p))))

;; should sort and do something about the leading "+"
(define (show-poly px)
  (foldr string-append "" 
    (foldr append '()   
      (map 
        (lambda (x) 
          (list
            " + "
            (number->string (cadr x)) 
            "x^" 
            (number->string (car x)))) 
        px))))

(define quadratic
  (lambda (a b c)
    (list (/ (- (* -1 b) (sqrt (- (* b b) (* 4 a c)))) (* 2 a))
          (/ (+ (* -1 b) (sqrt (- (* b b) (* 4 a c)))) (* 2 a)))))

(define (det2x2 a b c d)
  (poly-sum 
    (poly-prod (list '(1 -1) (list 0 a))
               (list '(1 -1) (list 0 d)))
    (list (list 0 (* -1 c b)))))


(define (root2x2 a b c d)
  (let* ([p (reverse (det2x2 a b c d))])
    (apply quadratic (map term-coef p))))

(define (to-poly xs c-ind)
  (define (row-iter cs indx)
    (cond 
      [(empty? cs) '()]
      [(= c-ind indx)
       (cons
         (list '(1  -1) (list 0 (car cs)))
         (row-iter (cdr cs) (add1 indx)))]
      [else 
        (cons 
          (list (list 0 (car cs)))
          (row-iter (cdr cs) (add1 indx)))]
      ))
  (row-iter xs 0))
 
(define (M-rI M)
  (map (lambda (ls n) (to-poly ls n)) M (range (length M))))



                  
(define (not-nth-term ls n)
  (if (= (length ls) (add1 n))
    (take ls n)
    (let ([hd (take ls  n)]
          [tl (take-right ls (- (length ls)  (add1 n)))])
      (append hd tl))))

(define (detpart r1 rns n)
  (let* ([p (list-ref r1 n)]
         [m (map (lambda (x) (not-nth-term x n)) rns)]
         [scalar (expt -1 n)]
         )
    (poly-prod (poly-scalar-prod scalar p) (det3x3 m))))
(define (det3x3 M)
  (let* ([len (length M)]
         [indices (range len)])
    (if (= (length M) 2)
      (poly-diff 
        (poly-prod (caar M) (cadadr M))
        (poly-prod (caadr M) (cadar M)))
      (let* ([row-1 (car M)]
             [rows-rest (cdr M)]
             )
        (foldr poly-sum '((0 0)) 
               (map (lambda (x) (detpart row-1 rows-rest x)) indices))))))
      
(define m '( (1 2)(4 5)) )
(define m1 '((((1 . -1) (0 . 1)) ((0 . 2))) (((0 . 4)) ((1 . -1) (0 . 5)))))
(define m2 '((((1    1) (0   2)) ((1   1))) (((1   1)) ((1    1)        ))))
(define m3 '(((1 1 2) (0 2 2) (-1 1 3)))) ;works!
(define m4 '(
             ((1 . 2 . 3)
              (4 . 5 . 6)
              (7 . 8 . 9))))
