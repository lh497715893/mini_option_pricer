#! /usr/bin/python
# -*- coding: utf8 -*-
from  Tkinter import *
import tkMessageBox
import ttk
import test

class option_pricer(Frame):

    def createWidget(self):

        content = Frame(self)
        input_lab = Label(content, text='Input', font=("Helvetica", 13))
        result = Label(content, text="Output:")
        # q1_result = Entry(content)######################

        content.grid(column=0, row=0)
        input_lab.grid(column=0, row=0, columnspan=5)
        result.grid(column=1, row=16)
        # q1_result.grid(column=1,row=18)################
        # q1_result.configure(state="disabled")################################
        Label(content, text=" ").grid(row=0, pady=20)
        self.resultContent = Label(content, text='')
        self.resultContent.grid(column=1, row=16, columnspan=3, pady=5)
        self.q1_parameter_p1 = Label(content, text='Spot Price of asset S1(S1)')
        self.q1_parameter_p1_layout = Entry(content)
        self.q1_parameter_p2 = Label(content, text='Spot Price of asset S2(S2)')
        self.q1_parameter_p2_layout = Entry(content)
        self.q2_parameter_p1 = Label(content, text='Volatility of asset S1(v1)')
        self.q2_parameter_p1_layout = Entry(content)
        self.q2_parameter_p2 = Label(content, text='Volatility of asset S2(v2)')
        self.q2_parameter_p2_layout = Entry(content)
        self.q3_parameter_p1 = Label(content, text='Time to maturity(T)')
        self.q3_parameter_p1_layout = Entry(content)
        self.q3_parameter_p2 = Label(content, text='Risk-free rate(r)')
        self.q3_parameter_p2_layout = Entry(content)
        self.q4_parameter_p1 = Label(content, text='Strike price(K)')
        self.q4_parameter_p1_layout = Entry(content)
        self.q4_parameter_p2 = Label(content, text='Option Type')
        self.optionType = StringVar()
        self.typeBtn1 = Radiobutton(content, text="CALL", variable=self.optionType, value="C")
        self.typeBtn2 = Radiobutton(content, text="PUT", variable=self.optionType, value="P")
        self.q5_parameter_p1 = Label(content, text='Correlation')
        self.q5_parameter_p1_layout = Entry(content)
        self.q6_parameter_p1 = Label(content, text='step No.')
        self.q6_parameter_p1_layout = Entry(content)
        self.q7_parameter_p1 = Label(content, text='path No.')
        self.q7_parameter_p1_layout = Entry(content)
        self.q7_parameter_p2 = Label(content, text='Control Variate')
        self.cvType = StringVar()
        self.param72cvBtn1 = Radiobutton(content, text="Std", variable=self.cvType, value="NULL")
        self.param72cvBtn2 = Radiobutton(content, text="Control Var", variable=self.cvType, value="CV")
        self.q8_parameter_p1 = Label(content, text='Repo rate')
        self.q8_parameter_p1_layout = Entry(content)
        self.q8_parameter_p2 = Label(content, text='Option premium')
        self.q8_parameter_p2_layout = Entry(content)
        self.showWidgets()
        self.typeBtn1.invoke()  # select Call Option as default
        self.param72cvBtn1.invoke()  # select NULL as defalt cv value

        Label(content, text=" ").grid(row=13, pady=10)

        self.questionValue = StringVar()


        self.players = ttk.Combobox(content, textvariable=self.questionValue)
        self.players["values"] = ("European call/put option", "Implied volatility calculator",
                                  "American call/put option", "GeoMetric Asian option",
                                  "Arithmetic Asian option",
                                  "GeoMetric basket option", "Arithmetic basket option")
        self.players.current(0)

        self.players.bind("<<ComboboxSelected>>", self.selected_info)
        self.players.grid(column=1, row=14,ipadx=18)

        calcBtn = Button(content, text='Result', command=self.calculate)
        calcBtn.grid(column=2, row=14)

    def selected_info(self,*args):
        a = self.questionValue.get()
        self.selected()

    def calculate(self):
        try:
            # Get Parameters
            selection = self.questionValue.get()
            S1 = float(self.q1_parameter_p1_layout.get())
            T = float(self.q3_parameter_p1_layout.get())
            r = float(self.q3_parameter_p2_layout.get())
            K = float(self.q4_parameter_p1_layout.get())
            type = self.optionType.get()
            # corr = float(self.param51.get())
            # n = int(self.param61.get())
            # path = int(self.param71.get())
            # cv = self.param72.get()
            # repo = float(self.param81.get())
            # trueValue = float(self.param82.get())
            # TODO: repo rate
            # repo = 0.3
            # trueValue = 10.0

            print 'TO DO-------execute calculation'
            resultPrice = 0.0000
            if selection == 'European call/put option':
                # repo = float(self.param81.get())
                sigma1 = float(self.q2_parameter_p1_layout.get())
                repo = 0.0
                if type == "C":
                    # stock, strike, time, maturity, volatility, repo, rfr
                    # resultPrice = optFunc.call_black_scholes(S1, K, 0.0, T, sigma1, repo, r
                    resultPrice = test.Block_Scholes_call(S1,K,T,sigma1,repo,r)
                elif type == "P":
                    # resultPrice = optFunc.put_black_scholes(S1, K, 0.0, T, sigma1, repo, r)
                    resultPrice = test.Block_Scholes_put(S1,K,T,sigma1,repo,r)
                else:
                    resultPrice = 0.0
                    # Q1
            elif selection == 'Implied volatility calculator':
                repo = float(self.q8_parameter_p1_layout.get())
                premium = float(self.q8_parameter_p2_layout.get())
                if type == 'C':
                    # S, K, t, T, q, r, pTrue
                    # resultPrice = optFunc.implied_vol_C(S1, K, 0.0, T, repo, r, trueValue)
                    resultPrice = test.vol_call(S1,K,0.2,T,r,premium)
                elif type == 'P':
                    # resultPrice = optFunc.implied_vol_P(S1, K, 0.0, T, repo, r, trueValue)
                    resultPrice = test.vol_put(S1,K,repo,T,r,premium)
                else:
                    resultPrice = 0.0
                    # Q2
            elif selection == "American call/put option":
                n = int(self.q6_parameter_p1_layout.get())
                sigma1 = float(self.q2_parameter_p1_layout.get())
                # S, K, r, T, sigma, N, type
                # resultPrice = optFunc.bino_tree(S1, K, r, T, sigma1, n, type)
                resultPrice = test.binoTreeMethod(S1,sigma1,r,T,K,n,type)
                #  Q3
            elif selection == 'GeoMetric Asian option':
                n = int(self.q6_parameter_p1_layout.get())
                sigma1 = float(self.q2_parameter_p1_layout.get())
                # S, sigma, r, t, K, n, type
                resultPrice = test.geoAsianOption(S1,sigma1,r,T,K,n,type)
                # Q4
            elif selection == 'Arithmetic Asian option':
                n = int(self.q6_parameter_p1_layout.get())
                path = int(self.q7_parameter_p1_layout.get())
                sigma1 = float(self.q2_parameter_p1_layout.get())
                cv = self.cvType.get()
                # S, sigma, r, T, K, step, type, path, cv
                resultPrice = test.arith_asian_option(S1,sigma1,r,T,K,n,type,path,cv)
                # Q5
            elif selection == 'GeoMetric basket option':
                S2 = float(self.q1_parameter_p2_layout.get())
                sigma1 = float(self.q2_parameter_p1_layout.get())
                sigma2 = float(self.q2_parameter_p2_layout.get())
                corr = float(self.q5_parameter_p1_layout.get())
                # S1, S2, sigma1, sigma2, r, T, K ,corr, type
                resultPrice = test.geoBasketOption(S1,S2,sigma1,sigma2,r,T,K,corr,type)
                # Q6
            elif selection == 'Arithmetic basket option':
                S2 = float(self.q1_parameter_p2_layout.get())
                sigma1 = float(self.q2_parameter_p1_layout.get())
                sigma2 = float(self.q2_parameter_p2_layout.get())
                corr = float(self.q5_parameter_p1_layout.get())
                path = int(self.q7_parameter_p1_layout.get())
                cv = self.cvType.get()
                # S1, S2, sigma1, sigma2, r, T, K, corr, type, path, cv
                # resultPrice = optFunc.arith_basket(S1, S2, sigma1, sigma2, r, T, K, corr, type, path, cv)
                resultPrice = test.arith_basket(S1,S2,sigma1,sigma2,r,T,K,corr,type,path,cv)
                # Q7
            # -------------------------------
            # show result in result box
            self.resultContent['text'] = "%.5f" % resultPrice
        except ValueError as e:
            print e
            tkMessageBox.showinfo("Please check your parameters", "Please check and input all parameters")

    def selected(self):
        # selection = self.questionValue.get()
        select = self.questionValue.get()
        if select=="":
            select = "European call/put option"
        # selection =
        # print "self.questionValue.get():"+str(selection)
        # print "qusetion", selection
        if select == 'European call/put option':
            self.showWidgets()
            self.q1_parameter_p1_layout.configure(state='normal')
            self.q1_parameter_p2_layout.configure(state='disabled')
            self.q2_parameter_p1_layout.configure(state='normal')
            self.q2_parameter_p2_layout.configure(state='disabled')
            self.q3_parameter_p1_layout.configure(state='normal')
            self.q3_parameter_p2_layout.configure(state='normal')
            self.q4_parameter_p1_layout.configure(state='normal')
            self.q5_parameter_p1_layout.configure(state='disabled')
            self.q6_parameter_p1_layout.configure(state='disabled')
            self.q7_parameter_p1_layout.configure(state='disabled')
            self.param72cvBtn1.configure(state='disabled')
            self.param72cvBtn2.configure(state='disabled')
            self.q8_parameter_p1_layout.configure(state='disabled')
            self.q8_parameter_p2_layout.configure(state='disabled')
            self.typeBtn1.configure(state='normal')
            self.typeBtn2.configure(state='normal')

            # Q1
        elif select == 'Implied volatility calculator':
            self.showWidgets()
            self.q1_parameter_p1_layout.configure(state='normal')
            self.q1_parameter_p2_layout.configure(state='disabled')
            self.q2_parameter_p1_layout.configure(state='disabled')
            self.q2_parameter_p2_layout.configure(state='disabled')
            self.q3_parameter_p1_layout.configure(state='normal')
            self.q3_parameter_p2_layout.configure(state='normal')
            self.q4_parameter_p1_layout.configure(state='normal')
            self.q5_parameter_p1_layout.configure(state='disabled')
            self.q6_parameter_p1_layout.configure(state='disabled')
            self.q7_parameter_p1_layout.configure(state='disabled')
            self.param72cvBtn1.configure(state='disabled')
            self.param72cvBtn2.configure(state='disabled')
            self.q8_parameter_p1_layout.configure(state='normal')
            self.q8_parameter_p2_layout.configure(state='normal')
            self.typeBtn1.configure(state='normal')
            self.typeBtn2.configure(state='normal')

        elif select == 'American call/put option':
            self.showWidgets()
            self.q1_parameter_p1_layout.configure(state='normal')
            self.q1_parameter_p2_layout.configure(state='disabled')
            self.q2_parameter_p1_layout.configure(state='normal')
            self.q2_parameter_p2_layout.configure(state='disabled')
            self.q3_parameter_p1_layout.configure(state='normal')
            self.q3_parameter_p2_layout.configure(state='normal')
            self.q4_parameter_p1_layout.configure(state='normal')
            self.q5_parameter_p1_layout.configure(state='disabled')
            self.q6_parameter_p1_layout.configure(state='normal')
            self.q7_parameter_p1_layout.configure(state='disabled')
            self.param72cvBtn1.configure(state='disabled')
            self.param72cvBtn2.configure(state='disabled')
            self.q8_parameter_p1_layout.configure(state='disabled')
            self.q8_parameter_p2_layout.configure(state='disabled')
            self.typeBtn1.configure(state='normal')
            self.typeBtn2.configure(state='normal')
            # Q3
        elif select == 'GeoMetric Asian option':
            self.showWidgets()
            self.q1_parameter_p1_layout.configure(state='normal')
            self.q1_parameter_p2_layout.configure(state='disabled')
            self.q2_parameter_p1_layout.configure(state='normal')
            self.q2_parameter_p2_layout.configure(state='disabled')
            self.q3_parameter_p1_layout.configure(state='normal')
            self.q3_parameter_p2_layout.configure(state='normal')
            self.q4_parameter_p1_layout.configure(state='normal')
            self.q5_parameter_p1_layout.configure(state='disabled')
            self.q6_parameter_p1_layout.configure(state='normal')
            self.q7_parameter_p1_layout.configure(state='disabled')
            self.param72cvBtn1.configure(state='disabled')
            self.param72cvBtn2.configure(state='disabled')
            self.q8_parameter_p1_layout.configure(state='disabled')
            self.q8_parameter_p2_layout.configure(state='disabled')
            self.typeBtn1.configure(state='normal')
            self.typeBtn2.configure(state='normal')
            # Q4
        elif select == 'Arithmetic Asian option':
            self.showWidgets()
            self.q1_parameter_p1_layout.configure(state='normal')
            self.q1_parameter_p2_layout.configure(state='disabled')
            self.q2_parameter_p1_layout.configure(state='normal')
            self.q2_parameter_p2_layout.configure(state='disabled')
            self.q3_parameter_p1_layout.configure(state='normal')
            self.q3_parameter_p2_layout.configure(state='normal')
            self.q4_parameter_p1_layout.configure(state='normal')
            self.q5_parameter_p1_layout.configure(state='disabled')
            self.q6_parameter_p1_layout.configure(state='normal')
            self.q7_parameter_p1_layout.configure(state='normal')
            self.param72cvBtn1.configure(state='normal')
            self.param72cvBtn2.configure(state='normal')
            self.q8_parameter_p1_layout.configure(state='disabled')
            self.q8_parameter_p2_layout.configure(state='disabled')
            self.typeBtn1.configure(state='normal')
            self.typeBtn2.configure(state='normal')
            # Q5
        elif select == 'GeoMetric basket option':
            self.showWidgets()
            self.q1_parameter_p1_layout.configure(state='normal')
            self.q1_parameter_p2_layout.configure(state='normal')
            self.q2_parameter_p1_layout.configure(state='normal')
            self.q2_parameter_p2_layout.configure(state='normal')
            self.q3_parameter_p1_layout.configure(state='normal')
            self.q3_parameter_p2_layout.configure(state='normal')
            self.q4_parameter_p1_layout.configure(state='normal')
            self.q5_parameter_p1_layout.configure(state='normal')
            self.q6_parameter_p1_layout.configure(state='disabled')
            self.q7_parameter_p1_layout.configure(state='disabled')
            self.param72cvBtn1.configure(state='disabled')
            self.param72cvBtn2.configure(state='disabled')
            self.q8_parameter_p1_layout.configure(state='disabled')
            self.q8_parameter_p2_layout.configure(state='disabled')
            self.typeBtn1.configure(state='normal')
            self.typeBtn2.configure(state='normal')
            # Q6
        elif select == 'Arithmetic basket option':
            # S1, S2, sigma1, sigma2, r, T, K, corr, type, path, cv
            self.showWidgets()
            self.q1_parameter_p1_layout.configure(state='normal')
            self.q1_parameter_p2_layout.configure(state='normal')
            self.q2_parameter_p1_layout.configure(state='normal')
            self.q2_parameter_p2_layout.configure(state='normal')
            self.q3_parameter_p1_layout.configure(state='normal')
            self.q3_parameter_p2_layout.configure(state='normal')
            self.q4_parameter_p1_layout.configure(state='normal')
            self.q5_parameter_p1_layout.configure(state='normal')
            self.q6_parameter_p1_layout.configure(state='disabled')
            self.q7_parameter_p1_layout.configure(state='normal')
            self.param72cvBtn1.configure(state='normal')
            self.param72cvBtn2.configure(state='normal')
            self.q8_parameter_p1_layout.configure(state='disabled')
            self.q8_parameter_p2_layout.configure(state='disabled')
            self.typeBtn1.configure(state='normal')
            self.typeBtn2.configure(state='normal')
            # Q7

    def showWidgets(self):
        self.q1_parameter_p1.grid(column=0, row=4)
        self.q1_parameter_p1_layout.grid(column=1, row=4)

        self.q1_parameter_p2.grid(column=0, row=5)
        self.q1_parameter_p2_layout.grid(column=1, row=5)

        self.q2_parameter_p1.grid(column=2, row=4)
        self.q2_parameter_p1_layout.grid(column=3, row=4)

        self.q2_parameter_p2.grid(column=2, row=5)
        self.q2_parameter_p2_layout.grid(column=3, row=5)

        self.q3_parameter_p1.grid(column=2, row=6)
        self.q3_parameter_p1_layout.grid(column=3, row=6)

        self.q3_parameter_p2.grid(column=0, row=7)
        self.q3_parameter_p2_layout.grid(column=1, row=7)

        self.q4_parameter_p1.grid(column=0, row=6)
        self.q4_parameter_p1_layout.grid(column=1, row=6)

        self.q4_parameter_p2.grid(column=2, row=7)
        self.typeBtn1.grid(column=3, row=7, sticky=W)
        self.typeBtn2.grid(column=3, row=7, sticky=E)

        self.q5_parameter_p1.grid(column=0, row=10)
        self.q5_parameter_p1_layout.grid(column=1, row=10)

        self.q6_parameter_p1.grid(column=0, row=9)
        self.q6_parameter_p1_layout.grid(column=1, row=9)

        self.q7_parameter_p1.grid(column=2, row=9)
        self.q7_parameter_p1_layout.grid(column=3, row=9)

        self.q7_parameter_p2.grid(column=2, row=8)
        self.param72cvBtn1.grid(column=3, row=8, sticky=W)
        self.param72cvBtn2.grid(column=3, row=8, sticky=E)

        self.q8_parameter_p1.grid(column=0, row=8)
        self.q8_parameter_p1_layout.grid(column=1, row=8)

        self.q8_parameter_p2.grid(column=2, row=10)
        self.q8_parameter_p2_layout.grid(column=3, row=10)

    # def createWidget(self):
    def __init__(self, master):
        root.minsize(width=700, height=500)
        Frame.__init__(self, master)
        self.pack()
        self.createWidget()
        self.selected()

# v = StringVar(root)
# v.set("Q1 European call/put option")
# om = OptionMenu(root, v, "Q2 Implied volatility calculator",
#                     "Q3 American call/put option", "Q4 GeoMetric Asian option",
#                     "Q5 Arithmetic Asian option", "Q6 GeoMetric basket option",
#                     "Q7 Arithmetic basket option")
# om.pack()
# om.grid(column=0,row=14)

root = Tk()
root.title("COMP7405 Assign GroupAL")
name = StringVar()
optFunc = option_pricer(master=root)
optFunc.mainloop()
