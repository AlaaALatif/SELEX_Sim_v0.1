import time
import random
import math
import numpy as np
from sklearn import linear_model
from sklearn.svm import SVR
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import explained_variance_score, mean_absolute_error, mean_squared_error, median_absolute_error, r2_score

##
sampleNum = 100

# generate fake data set for testing multiple mutlivariate regression
X = np.zeros((sampleNum, 2))
Y = np.zeros((sampleNum, 4))

for i, x in enumerate(X):
    X[i][0] = np.random.uniform(0, 1)
    X[i][1] = np.random.uniform(5, 10)
    Y[i][0] = (X[i][0]*2)+(X[i][1]*5)
    Y[i][1] = (X[i][0]*4)+(X[i][1]*7)
    Y[i][2] = X[i][0]+X[i][1]-5
    Y[i][3] = (X[i][0]*2)+(X[i][1]*2)-(X[i][0]*X[i][1]*2) # note that this target component is a non-linear function of x
    
# create training set
X_train = X[int(sampleNum*0.5):]
Y_train = Y[int(sampleNum*0.5):]

# create test set
X_test = X[:int(sampleNum*0.5)]
Y_test = Y[:int(sampleNum*0.5)]

# Machine learning: regression models
class Predictors:

    def decisionTreeRegression(self, x, y, xTest, yTest):
        model = DecisionTreeRegressor(max_depth=6)
        model.fit(x, y)
        preds = model.predict(xTest)
        print('Decision Tree Regression')
       # print('Coefficient: \n', model.coef_)
    # overall mean square error
        print("Residual sum-of-squares: %.2f" % np.mean((preds - yTest)**2))
    # overall max square error
        print("max residual sum-of-squares: %.2f" % np.amax((preds - yTest)**2))
    # overall variance score
        print("Variance score: %.2f" % explained_variance_score(yTest, preds, multioutput='uniform_average'))
    # Mean absolute error
        print("Mean absolute error: %.2f" % mean_absolute_error(yTest, preds))
    # Mean squared error
        print("Mean squared error: %.2f" % mean_squared_error(yTest, preds))
    # Median absolute error #THIS DOES NOT SUPPORT MULTIPLE OUTPUTS
        #print("Median absolute error: %.2f" % median_absolute_error(yTest, preds))
    # R^2 score
        print("R squared score: %.2f" % r2_score(yTest, preds))
        return model, yTest, preds





    
    def supportVectorRegression(self, x, y, xTest, yTest):

        #model = linear_model.BayesianRidge()
        models = [None for i in range(y.shape[1])]
        preds = np.zeros((yTest.shape[0], yTest.shape[1]))
        

        for i in range(y.shape[1]):
            models[i] = SVR()

            models[i].fit(x, y[:, i])

            preds[:, i] = models[i].predict(xTest)
        print('Support Vector Regression')
    # model learning rate
        #print(["learning rate: alpha = %.2f " % models[j].alpha_ for j in range(len(models))])
    # model coefficients
       # print('Coefficient: \n', model.coef_)
    # overall mean square error
        print("Residual sum-of-squares: %.2f" % np.mean((preds - yTest)**2))
    # overall max square error
        print("max residual sum-of-squares: %.2f" % np.amax((preds - yTest)**2))
    # overall variance score
       # print("Variance score: %.2f" % model.score(xTest, yTest))

        return models



## LINEAR MODELS

# Ordinary Least Squares
    def linearRegression(self, x, y, xTest, yTest):
    # object instantiation
        model = linear_model.LinearRegression()
    # model training
        model.fit(x, y)
    # predictions on test set
        preds = model.predict(xTest)

        print('Linear Regression')
    # model coefficients
        print('Coefficient: \n', model.coef_)
    # overall mean square error
        print("Residual sum-of-squares: %.2f" % np.mean((preds - yTest)**2))
    # overall max square error
        print("max residual sum-of-squares: %.2f" % np.amax((preds - yTest)**2))
    # overall variance score
        print("Variance score: %.2f" % model.score(xTest, yTest))

        for j in range(preds.shape[1]):
    # mean square error of first target component
            print("Residual sum-of-squares y0: %.2f" % np.mean((preds[:,j] - yTest[:,j])**2))
    # max square error of first target component
            print("max residual sum-of-squares y0: %.2f" % np.amax((preds[:,j] - yTest[:,j])**2))
    # variance score
        #print("Variance score: %.2f" % model.score(xTest, yTest[:,0]))
    # mean square error of first target component
        #print("Residual sum-of-squares y1: %.2f" % np.mean((preds[:,1] - yTest[:,1])**2))
    # max square error of second target component
        #print("max residual sum-of-squares y1: %.2f" % np.amax((preds[:,1] - yTest[:,1])**2))
    # variance score
        #print("Variance score: %.2f" % model.score(xTest, yTest[:,1]))
    # mean square error of third target component
        #print("Residual sum-of-squares y2: %.2f" % np.mean((preds[:,2] - yTest[:,2])**2))
    # max square error of third target component
        #print("max residual sum-of-squares y2: %.2f" % np.amax((preds[:,2] - yTest[:,2])**2))
    # variance score
        #print("Variance score: %.2f" % model.score(xTest, yTest[:,2]))
    # mean square error of fourth target component
        #print("Residual sum-of-squares y3: %.2f" % np.mean((preds[:,3] - yTest[:,3])**2))
    # max square error of fourth target component
        #print("max residual sum-of-squares y3: %.2f" % np.amax((preds[:,3] - yTest[:,3])**2))
    #  variance score
        #print("Variance score: %.2f" % model.score(xTest, yTest[:,3]))
        print j 
        return model

    def ridgeRegression(self, x, y, xTest, yTest):

        model = linear_model.RidgeCV(alphas = [0.1, 1.0, 10])

        model.fit(x, y)

        preds = model.predict(xTest)
        
        print('Ridge Regression')

    # model learning rate
        print("learning rate: alpha = %.2f " % model.alpha_ )
    # model coefficients
        print('Coefficient: \n', model.coef_)
    # overall mean square error
        print("Residual sum-of-squares: %.2f" % np.mean((preds - yTest)**2))
    # overall max square error
        print("max residual sum-of-squares: %.2f" % np.amax((preds - yTest)**2))
    # overall variance score
        print("Variance score: %.2f" % model.score(xTest, yTest))


    # mean square error of first target component
        print("Residual sum-of-squares y0: %.2f" % np.mean((preds[:,0] - yTest[:,0])**2))
    # max square error of first target component
        print("max residual sum-of-squares y0: %.2f" % np.amax((preds[:,0] - yTest[:,0])**2))
    # variance score
        #print("Variance score: %.2f" % model.score(xTest, yTest[:,0]))
    # mean square error of first target component
        print("Residual sum-of-squares y1: %.2f" % np.mean((preds[:,1] - yTest[:,1])**2))
    # max square error of second target component
        print("max residual sum-of-squares y1: %.2f" % np.amax((preds[:,1] - yTest[:,1])**2))
    # variance score
        #print("Variance score: %.2f" % model.score(xTest, yTest[:,1]))
    # mean square error of third target component
        print("Residual sum-of-squares y2: %.2f" % np.mean((preds[:,2] - yTest[:,2])**2))
    # max square error of third target component
        print("max residual sum-of-squares y2: %.2f" % np.amax((preds[:,2] - yTest[:,2])**2))
    # variance score
        #print("Variance score: %.2f" % model.score(xTest, yTest[:,2]))
    # mean square error of fourth target component
        print("Residual sum-of-squares y3: %.2f" % np.mean((preds[:,3] - yTest[:,3])**2))
    # max square error of fourth target component
        print("max residual sum-of-squares y3: %.2f" % np.amax((preds[:,3] - yTest[:,3])**2))
    #  variance score
        #print("Variance score: %.2f" % model.score(xTest, yTest[:,3]))
    
        return model




    def bayesianRidgeRegression(self, x, y, xTest, yTest):

        #model = linear_model.BayesianRidge()
        models = [None for i in range(y.shape[1])]
        preds = np.zeros((yTest.shape[0], yTest.shape[1]))
        

        for i in range(y.shape[1]):
            models[i] = linear_model.BayesianRidge()

            models[i].fit(x, y[:, i])

            preds[:, i] = models[i].predict(xTest)
        print('Bayesian Ridge Regression')
    # model learning rate
        print(["learning rate: alpha = %.2f " % models[j].alpha_ for j in range(len(models))])
    # model coefficients
       # print('Coefficient: \n', model.coef_)
    # overall mean square error
        print("Residual sum-of-squares: %.2f" % np.mean((preds - yTest)**2))
    # overall max square error
        print("max residual sum-of-squares: %.2f" % np.amax((preds - yTest)**2))
    # overall variance score
       # print("Variance score: %.2f" % model.score(xTest, yTest))


    # mean square error of first target component
        print("Residual sum-of-squares y0: %.2f" % np.mean((preds[:,0] - yTest[:,0])**2))
    # max square error of first target component
        print("max residual sum-of-squares y0: %.2f" % np.amax((preds[:,0] - yTest[:,0])**2))
    # variance score
        #print("Variance score: %.2f" % model.score(xTest, yTest[:,0]))
    # mean square error of first target component
        print("Residual sum-of-squares y1: %.2f" % np.mean((preds[:,1] - yTest[:,1])**2))
    # max square error of second target component
        print("max residual sum-of-squares y1: %.2f" % np.amax((preds[:,1] - yTest[:,1])**2))
    # variance score
        #print("Variance score: %.2f" % model.score(xTest, yTest[:,1]))
    # mean square error of third target component
        print("Residual sum-of-squares y2: %.2f" % np.mean((preds[:,2] - yTest[:,2])**2))
    # max square error of third target component
        print("max residual sum-of-squares y2: %.2f" % np.amax((preds[:,2] - yTest[:,2])**2))
    # variance score
        #print("Variance score: %.2f" % model.score(xTest, yTest[:,2]))
    # mean square error of fourth target component
        print("Residual sum-of-squares y3: %.2f" % np.mean((preds[:,3] - yTest[:,3])**2))
    # max square error of fourth target component
        print("max residual sum-of-squares y3: %.2f" % np.amax((preds[:,3] - yTest[:,3])**2))
    #  variance score
        #print("Variance score: %.2f" % model.score(xTest, yTest[:,3]))
    
        return models



#TEST AREA
#regPred = Predictors()

#olsModel = regPred.linearRegression(X_train, Y_train, X_test, Y_test)
#ridgeModel = regPred.ridgeRegression(X_train, Y_train, X_test, Y_test)
#bayesModel = regPred.bayesianRidgeRegression(X_train, Y_train, X_test, Y_test)



