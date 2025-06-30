Google colab link: https://colab.research.google.com/drive/1ZFgGIhVFuunZtIllH1kUleCFePZYWsD8#scrollTo=1wiLUQ24PeUA

!pip install joblib
!pip install -U scikit-learn

# Import libraries
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.metrics import r2_score, mean_squared_error
import joblib
import matplotlib.pyplot as plt
import seaborn as sns

# Upload your training CSV file
from google.colab import files
uploaded = files.upload()

# Load training data
training_file_path = list(uploaded.keys())[0]
training_data = pd.read_csv(training_file_path)
print("✅ Training data loaded successfully!")

# Define input features and target
features = [
    "nucleotide_composition_A", "nucleotide_composition_U",
    "nucleotide_composition_G", "nucleotide_composition_C",
    "gc_content", "minimum_free_energy", "Molecular Weight (g/mol)",
    "Hydrogen Bond Donors", "Hydrogen Bond Acceptors", "LogP", "TPSA (Å²)",
    "Atom Count", "Chiral Atom Count", "num_electrostatic_contacts",
    "avg_electrostatic_distance", "num_vdw_contacts", "avg_vdw_distance"
]
target = "Binding Affinity (kcal/mol)"

# Split data
X_train = training_data[features]
y_train = training_data[target]

# Imputation and scaling
imputer = SimpleImputer(strategy='mean')
scaler = StandardScaler()
target_scaler = StandardScaler()

X_train_imputed = imputer.fit_transform(X_train)
X_train_scaled = scaler.fit_transform(X_train_imputed)
y_train_scaled = target_scaler.fit_transform(y_train.values.reshape(-1, 1)).ravel()

rf_model = RandomForestRegressor(random_state=42)
param_grid_rf = {
    'n_estimators': [100, 200],
    'max_depth': [10, 20, None],
    'min_samples_split': [5, 10],
    'min_samples_leaf': [2, 4],
    'max_features': ['sqrt', 'log2']
}

grid_search_rf = GridSearchCV(rf_model, param_grid_rf, cv=5, scoring='r2', n_jobs=-1)
grid_search_rf.fit(X_train_scaled, y_train_scaled)

best_rf_model = grid_search_rf.best_estimator_
print("✅ Model trained. Best parameters:", grid_search_rf.best_params_)

y_pred_train = best_rf_model.predict(X_train_scaled)
y_pred_train_unscaled = target_scaler.inverse_transform(y_pred_train.reshape(-1, 1)).flatten()
y_train_unscaled = target_scaler.inverse_transform(y_train_scaled.reshape(-1, 1)).flatten()

print("🔎 R2 Score:", r2_score(y_train_unscaled, y_pred_train_unscaled))
mse = mean_squared_error(y_train_unscaled, y_pred_train_unscaled)
rmse = np.sqrt(mse)
print("🔎 RMSE:", rmse)

pipeline = {
    'model': best_rf_model,
    'imputer': imputer,
    'scaler': scaler,
    'target_scaler': target_scaler,
    'features': features
}
joblib.dump(pipeline, 'binding_affinity_pipeline.pkl')
print("✅ Full pipeline saved as 'binding_affinity_pipeline.pkl'")

# Upload new PDB features CSV file
uploaded = files.upload()
new_pdb_file = list(uploaded.keys())[0]
new_data = pd.read_csv(new_pdb_file)
pipeline = joblib.load('binding_affinity_pipeline.pkl')

# Prediction function
def predict_binding_affinity(new_data, pipeline):
    try:
        required_features = pipeline['features']
        missing = [f for f in required_features if f not in new_data.columns]
        if missing:
            raise ValueError(f"Missing required features: {missing}")

        X_new = new_data[required_features]
        X_imputed = pipeline['imputer'].transform(X_new)
        X_scaled = pipeline['scaler'].transform(X_imputed)
        y_scaled_pred = pipeline['model'].predict(X_scaled)
        y_pred = pipeline['target_scaler'].inverse_transform(y_scaled_pred.reshape(-1, 1)).flatten()

        result_df = pd.DataFrame({
            "PDB ID": new_data.get("PDB ID", [f"PDB_{i}" for i in range(len(y_pred))]),
            "Predicted Binding Affinity": y_pred
        })

        output_filename = f"Predicted_Binding_Affinity_{new_pdb_file}"
        result_df.to_csv(output_filename, index=False)
        print(f"✅ Predictions saved to: {output_filename}")
        return result_df
    except Exception as e:
        print(f"❌ Error: {e}")

# Run prediction
predicted_df = predict_binding_affinity(new_data, pipeline)
predicted_df.head()

def plot_feature_importance(model, feature_names):
    importances = model.feature_importances_
    sorted_idx = np.argsort(importances)[::-1]
    plt.figure(figsize=(10, 6))
    sns.barplot(x=importances[sorted_idx], y=np.array(feature_names)[sorted_idx])
    plt.title("Random Forest Feature Importances")
    plt.tight_layout()
    plt.show()

plot_feature_importance(pipeline['model'], pipeline['features'])
