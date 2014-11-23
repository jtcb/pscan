(ns pscan.core 
  (:require [clojure.core.reducers :as r]
            [clojure.set :as s]))

(def pts [[0 0] [0 1] [0 2] [1 0] [1 1] [1 2] [2 0] [2 1] [2 2]])

(defn r2-dist
  "Euclidian distance for points in R^2 (for testing)
   pt is a vector of length 2"
  [a b]
  (if (or (nil? a) (nil? b))
    Double/POSITIVE_INFINITY
  ;else
    (let [dx (- (a 0) (b 0))
          dy (- (a 1) (b 1))]
      (Math/sqrt (+ (* dx dx) (* dy dy))))))

(defn region-query-serial
  "Serial region query; returns all pts within eps of pt by metric."
  [pt pts metric eps]
  (letfn [(is-local? [x] (<= (metric pt x) eps))]
    (into [] (filter is-local? pts))))


(comment
(def x (into [] (take 10000000 (repeatedly rand))))

(time (into [] (filter #(> % 0.999999) x)))

(time (into #{} (filter #(> % 0.999999) x)))

(time (r/fold (r/monoid into vector) conj (r/filter #(> % 0.999999) x)))

(time (r/fold (r/monoid into hash-set) conj (r/filter #(> % 0.999999) x)))

(time (r/fold (r/monoid s/union hash-set) conj (r/filter #(> % 0.999999) x)))

(region-query-serial [0 0] pts r2-dist 2)
  
  )

