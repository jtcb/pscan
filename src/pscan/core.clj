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

(defn region-query
  "Serial region query; returns all pts within eps of pt by metric."
  [pt pts metric eps]
  (letfn [(is-local? [x] (<= (metric pt x) eps))]
    (into [] (filter is-local? pts))))

(defn expand-cluster
  "Expand cluster at pt
   Returns false if there is not cluster to be had
   Returns the core points of the cluster and an updated visited set otherwise" 
  [pt pts metric eps min-pts visited]

  (let [local (region-query pt pts metric eps)]
    (if (< (count local) min-pts)
      false
    ;else
      (loop [s (disj (set (remove #(visited %) local)) pt)
             c #{pt} 
             v (conj visited pt)]
        (if (empty? s)
          {:core c, :visited v} 
        ;else 
          (let [current-point (first s)
                result (region-query current-point pts metric eps)
                size (count result)
                unvisited-result (filter #(not (v %)) result)]
                
            (if (>= size min-pts) ;new core point
              (recur (disj (into s unvisited-result) current-point)
                     (conj c current-point)
                     (conj v current-point))
            ;else
              (recur (disj s current-point)
                     c 
                     (conj v current-point)))))))))





(comment

(def star [[0 0] [0 1] [1 0] [-1 0] [0 -1]])

(expand-cluster [0 0] star r2-dist 1 2 #{})

(expand-cluster [0 0] star r2-dist 1 5 #{})

(expand-cluster [-1 0] star r2-dist 1.5 5 #{})

(def x (into [] (take 10000000 (repeatedly rand))))

(time (into [] (filter #(> % 0.999999) x)))

(time (into #{} (filter #(> % 0.999999) x)))

(time (r/fold (r/monoid into vector) conj (r/filter #(> % 0.999999) x)))

(time (r/fold (r/monoid into hash-set) conj (r/filter #(> % 0.999999) x)))

(time (r/fold (r/monoid s/union hash-set) conj (r/filter #(> % 0.999999) x)))

(region-query [0 0] pts r2-dist 2)
  
(map #(region-query % pts r2-dist 1) '([0 0] [0 1]))
)

