import shap
# model = xgboost.XGBRegressor().fit(X_train, y_train)

# explain the model's predictions using SHAP
# (same syntax works for LightGBM, CatBoost, scikit-learn, transformers, Spark, etc.)
def cal_shap_values(model, data):
    explainer = shap.TreeExplainer(model)
    shap_values = explainer(data)
    return explainer, shap_values




