(ns pscan.core 
  (:require [clojure.core.reducers :as r]
            [clojure.set :as s]
            [clojure.string :refer [split]]
            [clj-biosequence.core :as cbs]
            [clojure.java.io :as io])
  (:import match.Match2))

(defn r2-dist
  "NegativeEuclidian distance for points in R^2 (for testing)
   pt is a vector of length 2"
  [a b]
  (if (or (nil? a) (nil? b))
    Double/NEGATIVE_INFINITY
  ;else
    (let [dx (- (a 0) (b 0))
          dy (- (a 1) (b 1))]
      (-
        (Math/sqrt (+ (* dx dx) (* dy dy)))))))

(defn region-query-s
  "Serial region query; returns all pts *at least* eps similar to pt."
  [pt pts metric eps]
  (letfn [(is-local? [x] (>= (metric pt x) eps))]
    (into [] (filter is-local? pts))))

(defn region-query-p
  "Parallel region query; returns all pts at least eps similar to pt."
  [pt pts metric eps]
  (letfn [(is-local? [x] (>= (metric pt x) eps))]
    (r/fold 5 (r/monoid into vector) conj (r/filter is-local? pts))))

;; TODO; find optimal subdivision for r/fold; memoize

;region-query-m memoized

;region-query-mp memoized, parallel

;; alias the appropriate region-query variant
(def region-query region-query-p)

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
                unvisited-result (remove #(v %) result)]
                
            (if (>= size min-pts) ;new core point
              (recur (disj (into s unvisited-result) current-point)
                     (conj c current-point)
                     (conj v current-point))
            ;else
              (recur (disj s current-point)
                     c 
                     (conj v current-point)))))))))

(defn dbscan
  "Run dbscan"
  [pts metric eps min-pts]
  
  (loop [points pts
         visited #{}
         clusters #{}]
    (let [pt (first points)]
      (cond
        (empty? points)
          clusters
        (visited pt) 
          (recur (rest points)
                 (conj visited pt)
                 clusters)
        :else
          (let [{new-cluster :core, new-visited :visited}
                (expand-cluster pt pts metric eps min-pts visited)]
            (recur (rest points)
                   new-visited
                   (conj clusters new-cluster))))))) 

(defn parse-float
  "Parse string s as float"
  [s]
  (Float/parseFloat s))

(defn blosum-score
  "Returns the BLOSUM score of a aligned to b
   BLOSUM62 matrix, gap penalty 11, extension penalty 1 (BLAST defaults)"
  [a b]
  (match.Match2/blosum62 a b))

(defn process-fasta
  "Convert clj_biosequence.core.fastaSequence into a keyword->string map"
  [x file]
  {:sequence (apply str (:sequence x)),
   :file file,
   :description (:description x)}) 

(defn read-fasta
  "Read protein data into a vector of maps"
  [file]
  (with-open [r (cbs/bs-reader (cbs/init-fasta-file file :iupacAminoAcids))]
    (let [fastas (cbs/biosequence-seq r)]
      (into [] (map #(process-fasta % file) fastas)))))

(defn protein-metric
  "blosum score between 2 fastaSequences
   Pseudo-metric (can be negative); larger is closer"
  [a b]
  (blosum-score (:sequence a) (:sequence b)))

(defn memoize-metric
  "Return a memoized version of metric f (symmetric)"
  [f]
  (let [mem (atom {})]
    (fn [a b]
      (if-let [e (find @mem #{a b})]
        (val e)
        (let [ret (f a b)]
          (swap! mem assoc #{a b} ret)
          ret)))))

;;;

(def files
  ["UniRef90_L7WSQ0.fa",  "UniRef90_Q2FZQ3.fa", "UniRef90_A7X1N2.fa",
   "UniRef90_L7WXN8.fa",  "UniRef90_Q5HMF1.fa", "UniRef90_C6F1E2.fa",
   "UniRef90_P0CK84.fa",  "UniRef90_Q5HPT5.fa", "UniRef90_I3VED1.fa",
   "UniRef90_Q0A457.fa"])

(def prefixed-files
  (mapv #(str "resources/" %) files))

(def sample-proteins
  (reduce into 
          (map read-fasta (take 3 prefixed-files)))) 

(time (def trial (dbscan sample-proteins protein-metric 100 4)))


(comment

(read-fasta "resources/simple.fasta")

(def pts [[0 0] [0 1] [0 2] [1 0] [1 1] [1 2] [2 0] [2 1] [2 2]])

(def star [[0 0] [0 1] [1 0] [-1 0] [0 -1]])

(def blobs [[0 0] [0 1] [1 0] [1 1]
            [100 100] [101 100] [100 101] [101 101]])

(dbscan blobs r2-dist 2 3)

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

