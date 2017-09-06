protein_complex_maps
====================

Kdrew's scripts for handling protein complex map data

How to run:
1. Generate elution profiles  (preprocessing)
2. Calculate correlations (feature extraction)
3. Convert to pairwise files (feature extraction)
4. Build feature matrix (feature extraction)
4b. Ensure Common ID (feature extraction)
5. Split benchmark (Model Training/make_benchmark)
6. Train classifier (Model Training)
7. Predict interactions (Model Training)
8. Evaluate predicted interactions (Evaluation)
9. Cluster interactions (Clustering)
10. Evaluate predicted clusters (Evaluation)


Preprocessing  (src/preprocessing_util/)
Feature Extraction (src/features/)
Model Training (src/model_fitting/)
    SVM, LDA, tpot, other machine learning
Clustering (src/clustering/)
Evaluation (src/evalution/)


Filename conventions:




