# Project Five: An ML-Powered Stock Predictor

This was an exploratory project into the world of machine learning and financial markets. The objective was to build a model that could predict short-term stock price movements based on historical data and sentiment analysis from news articles.

## The Challenge

Financial markets are notoriously difficult to predict. The challenge was to build a model that could find meaningful patterns in noisy data and provide predictions that were better than random chance.

## Our Approach

- **Data Collection:** We collected historical stock price data using the Alpha Vantage API and scraped financial news articles from various sources.
- **Feature Engineering:** We created a variety of features, including technical indicators (e.g., moving averages, RSI) and sentiment scores from the news articles (using NLP).
- **Model Selection:** We experimented with several models, including LSTM (Long Short-Term Memory) networks and Gradient Boosting machines.
- **Backtesting:** We rigorously backtested our models on historical data to evaluate their performance.

## Results and Learnings

While the model did not discover a "holy grail" for stock market prediction, it did show a slight predictive edge over a random baseline. The most significant learning was the difficulty of working with financial data and the importance of avoiding overfitting.

The sentiment analysis component proved to be particularly interesting, showing a correlation between positive news sentiment and short-term price increases.

This project was a fascinating introduction to the application of machine learning in a complex, real-world domain. It underscored the importance of domain knowledge, careful feature engineering, and robust validation techniques.
