import sys
import pandas as pd
from src.eur_predictor import EURPredictor

if __name__ == '__main__':
    data_file = sys.argv[1] if len(sys.argv) > 1 else 'sample_well.csv'
    days = int([arg for arg in sys.argv if arg.startswith('--days=')][0].split('=')[1]) if any(arg.startswith('--days=') for arg in sys.argv) else 45
    df = pd.read_csv(data_file)
    predictor = EURPredictor()
    result = predictor.predict(df, early_days=days)
    result.save_outputs()
    print('Forecast complete! Check outputs/ folder for Excel and plots.')