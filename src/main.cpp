/**
 * @file main.cpp
 * @author Zach Silberstein (zach.silberstein@gmail.com)
 * @brief Dynamical simulation and control of the unstable double inverted
 * pendulum (DIP) on a cart.
 * @date 2024-05-31
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <iostream>
#include "raylib.h"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"
#include "Eigen/Dense"
#include "DIP.h"
#include "LQR.h"

/**
 * 
 * @struct GUISlider
 * @brief Stores data for a GUI slider.
 * @param label Label above the slider.
 * @param currVale The current value of the slider.
 * @param minValue The minumum value of the slider.
 * @param maxValue The maximum value of the slider.
 * 
 */
struct GUISlider {
    const char *label;
    float currValue;
    float minValue;
    float maxValue;
};

typedef enum
{
    OPEN_LOOP = 0,
    CLQR,
    DLQR
} controlMode;

int main(void) {

    // Set up window
    const int screenWidth = 1600;
    const int screenHeight = 900;
    InitWindow(screenWidth, screenHeight, "DIP-Control");

    const int FPS = 100;
    SetTargetFPS(FPS); // Time between points in sim is 1/FPS
    const double dt = 1.0 / FPS;

    // Set style
    GuiSetStyle(DEFAULT, TEXT_SIZE, 17);

    // Conversion factor from meters to pixels
    const float meter = std::min(screenWidth, screenHeight) * 0.1f;

    // Cart dimensions in meters
    const Vector2 cartDims = {meter, meter * 0.5f};
    const float wheelRadius = meter * 0.2f;
    // Center of cart body in pixels
    const Vector2 cartTrackCenter = {screenWidth * 0.75f / 2.0f - cartDims.x / 2.0f,
                                     screenHeight * 0.55f - cartDims.y - wheelRadius};

    // X, Y coordinates of DIP system in pixels
    Vector2 cartPos;
    Vector2 leftWheelPos;
    Vector2 rightWheelPos;
    Vector2 massOnePos;
    Vector2 massTwoPos;

    // Sliders for initialization (label, initial value, min value, max value)
    GUISlider initSliders[] = {{"Cart Mass: %.2f kgs", 3.0f, 1.0f, 10.0f},
                               {"Mass 1: %.2f kgs", 1.0f, 0.5f, 5.0f},
                               {"Mass 2: %.2f kgs", 1.0f, 0.5f, 5.0f},
                               {"Length 1: %.2f m", 1.0f, 0.1f, 2.5f},
                               {"Length 2: %.2f m", 1.0f, 0.1f, 2.5f},
                               {"Cart Damping: %.2f Ns/m", 0.0f, 0.0f, 5.0f},
                               {"Theta 1 Damping: %.2f Nms/rad", 0.0f, 0.0f, 3.0f},
                               {"Theta 2 Damping: %.2f Nms/rad", 0.0f, 0.0f, 3.0f},
                               {"Theta 1: %.2f rads", 0.0f, -M_PIf, M_PIf},
                               {"Theta 2: %.2f rads", 0.0f, -M_PIf, M_PIf},
                               {"Cart Velocity: %.2f m/s", 0.0f, -2.0f, 2.0f},
                               {"Link 1 Angular Velocity: %.2f rads/s", 0.0f, -1.0f, 1.0f},
                               {"Link 2 Angular Velocity: %.2f rads/s", 0.0f, -2.0f, 2.0f}};

    // Sliders to display the state (label, initial value, min value, max value)
    GUISlider stateSliders[] = {{"Cart Position: %.2f m", 0.0f, -screenWidth * 0.325 / meter, screenWidth * 0.325 / meter},
                                {"Theta 1: %.2f rads", 0.0f, -M_PIf, M_PIf},
                                {"Theta 2: %.2f rads", 0.0f, -M_PIf, M_PIf},
                                {"Cart Velocity: %.2f m/s", 0.0f, -15.0f, 15.0f},
                                {"Link 1 Angular Velocity: %.2f rads/s", 0.0f, -25.0f, 25.0f},
                                {"Link 2 Angular Velocity: %.2f rads/s", 0.0f, -25.0f, 25.0f}};

    // Number of sliders to show, changes with control mode
    size_t numSliders = 13;

    // Location of start button, changes with control mode
    float startBtnHeight = screenHeight * 0.85;
    
    // Start button value
    bool startBtn = false;

    // Pause/resume button value
    bool pauseBtn = false;
    // Text for pause/resume button
    const char *btnText;

    // Dropdown menu vars
    int controlModeActive = OPEN_LOOP;
    bool controlModeEdit = false;

    // Reset button value
    bool resetBtn = false;

    // DIP object to be created and used later
    DIP* dip = nullptr;

    // LQR object to be created and used later
    LQR* lqr = nullptr;

    // Time in sim
    double time{0};

    // Store the state of the system
    Eigen::Vector<double, 6> state;

    // Main loop consists of initialization loop and simulation loop
    while (!WindowShouldClose()) {

        // Initialization loop
        while (!startBtn && !WindowShouldClose()) {

            // Update postion of the DIP system
            cartPos = cartTrackCenter;  // While not started, cart stays at center
            leftWheelPos = {cartPos.x, cartPos.y + cartDims.y};
            rightWheelPos = {cartPos.x + cartDims.x, cartPos.y + cartDims.y};
            massOnePos = {cartPos.x + meter * initSliders[3].currValue * sinf(initSliders[8].currValue) + cartDims.x / 2.0f,
                        cartPos.y - meter * initSliders[3].currValue * cosf(initSliders[8].currValue) + cartDims.y / 2.0f};
            massTwoPos = {massOnePos.x + meter * initSliders[4].currValue * sinf(initSliders[8].currValue + initSliders[9].currValue),
                        massOnePos.y - meter * initSliders[4].currValue * cosf(initSliders[8].currValue +  initSliders[9].currValue)};

            BeginDrawing();

            ClearBackground(WHITE);

            // Cart Track
            DrawLineEx((Vector2){screenWidth * 0.05, screenHeight * 0.55},
                    (Vector2){screenWidth * 0.7, screenHeight * 0.55}, 
                    5,
                    Fade(BLACK, 0.8f));
            // Cart (Darker with increasing mass)
            DrawRectangleV(cartPos, cartDims,
                        Fade(BLUE, initSliders[0].currValue / 20.0f + 0.5f));
            // Left wheel
            DrawCircleV(leftWheelPos, wheelRadius, Fade(RED, 0.95f));
            // Right wheel
            DrawCircleV(rightWheelPos, wheelRadius, Fade(RED, 0.95f));
            // Link 1
            DrawLineEx((Vector2){cartPos.x + cartDims.x / 2, cartPos.y + cartDims.y / 2},
                        massOnePos, 10, Fade(BLACK, 0.9f));
            // Link 2
            DrawLineEx(massOnePos, massTwoPos, 10, Fade(BLACK, 0.9f));
            // Mass 1 (Radius grows with mass)
            DrawCircleV(massOnePos, (sqrt(initSliders[1].currValue) + 1.0) * meter / 10, RED);
            // Mass 2 (Radius grows with mass)
            DrawCircleV(massTwoPos, (sqrt(initSliders[2].currValue) + 1.0) * meter / 10, RED);

            // Make GUI area on right quarter of window
            DrawLine(screenWidth * 0.75, 0, screenWidth * 0.75, screenWidth,
                    Fade(LIGHTGRAY, 0.6f));
            DrawRectangle(screenWidth * 0.75, 0, screenWidth * 0.25, screenHeight,
                        Fade(LIGHTGRAY, 0.3f));

            // Lock GUI if editing dropdown menu
            if (controlModeEdit) {
                GuiLock();
            }

            // Add GUI sliders
            for (size_t i = 0; i < numSliders; ++i) {
                GuiLabel((Rectangle){screenWidth * 0.775,
                                    screenHeight * (0.07f + i * 0.06f), 1400, 24},
                                    TextFormat(initSliders[i].label, initSliders[i].currValue));
                GuiSliderBar((Rectangle){screenWidth * 0.775,
                                        screenHeight * (0.09f + i * 0.06f),
                                        screenWidth * 0.2, screenHeight * 0.025},
                            NULL, NULL, &initSliders[i].currValue, initSliders[i].minValue,
                            initSliders[i].maxValue);
            }

            // Add start button
            if (GuiButton((Rectangle){screenWidth * 0.775, startBtnHeight,
                                    screenWidth * 0.2, screenHeight * 0.025},
                        "Start")) {
                startBtn = true;
            }

            GuiUnlock();

            // Control mode Dropdown
            GuiLabel(
                (Rectangle){screenWidth * 0.775, screenHeight * 0.01, 140, 24},
                "Control Mode");
            if (GuiDropdownBox(
                    (Rectangle){screenWidth * 0.775, screenHeight * 0.03,
                                screenWidth * 0.25 * 0.8, screenHeight * 0.025},
                    "Open Loop;Continuos LQR; Discrete LQR", &controlModeActive,
                    controlModeEdit)) {
                        controlModeEdit = !controlModeEdit;
                        if (controlModeActive == OPEN_LOOP) {
                            initSliders[8].minValue = -M_PIf;
                            initSliders[9].minValue = -M_PIf;
                            initSliders[8].maxValue = M_PIf;
                            initSliders[9].maxValue = M_PIf;
                            numSliders = 13;
                            startBtnHeight = screenHeight * 0.85;
                        }
                        else {
                            // Limit angles to be within small angle
                            // approximation when in control mode
                            initSliders[8].minValue = -M_PIf / 20;
                            initSliders[9].minValue = -M_PIf / 20;
                            initSliders[8].maxValue = M_PIf / 20;
                            initSliders[9].maxValue = M_PIf / 20;
                            numSliders = 10;
                            startBtnHeight = screenHeight * 0.67f;
                        }
                    }

            EndDrawing();
        }

        // Create DIP object with inital values
        if (!WindowShouldClose()) {
            dip = new DIP(initSliders[0].currValue, initSliders[1].currValue,
                        initSliders[2].currValue, initSliders[3].currValue,
                        initSliders[4].currValue, initSliders[5].currValue,
                        initSliders[6].currValue, initSliders[7].currValue, dt,
                        Eigen::Vector<double, 6>{
                            0, initSliders[8].currValue, initSliders[9].currValue,
                            initSliders[10].currValue, initSliders[11].currValue,
                            initSliders[12].currValue});

            // Create LQR object, let Q and R be identity matrix for now
            if (controlModeActive == CLQR) {
                lqr = new LQR(dip->getAc(), dip->getBc(),
                              Eigen::MatrixXd::Identity(6, 6),
                              Eigen::MatrixXd::Identity(1, 1), true);
            }
            else if (controlModeActive == DLQR) {
                lqr = new LQR(dip->getAd(), dip->getBd(),
                              Eigen::MatrixXd::Identity(6, 6),
                              Eigen::MatrixXd::Identity(1, 1), false);
            }
        }

        // Simulation loop
        while (!resetBtn && !WindowShouldClose()) {

            // Do not update if paused
            if (!pauseBtn) {
                
                // Update state of DIP
                if (controlModeActive == OPEN_LOOP) {
                     state = dip->updateState(0);
                }
                else {
                    // Update state with feedback force based on current state
                    state = dip->updateState(lqr->getFeedbackControl(dip->getState())(0));
                }                

                // Update state sliders
                for (size_t i = 0; i < std::size(state); ++i) {
                    stateSliders[i].currValue = state[i];
                }

                // Update time
                time += dt;

                // Update postion of DIP in pixels
                cartPos = {cartTrackCenter.x + ((float)state[0]) * meter,
                        cartTrackCenter.y};
                leftWheelPos = {cartPos.x, cartPos.y + cartDims.y};
                rightWheelPos = {cartPos.x + cartDims.x, cartPos.y + cartDims.y};
                massOnePos = {cartTrackCenter.x + meter * ((float)dip->getMassOnePos()[0]) +
                                    cartDims.x / 2.0f,
                            cartTrackCenter.y - meter * ((float)dip->getMassOnePos()[1]) +
                                    cartDims.y / 2.0f};
                massTwoPos = {cartTrackCenter.x + meter * ((float)dip->getMassTwoPos()[0]) +
                                    cartDims.x / 2.0f,
                            cartTrackCenter.y - meter * ((float)dip->getMassTwoPos()[1]) +
                                    cartDims.y / 2.0f};

            }

            BeginDrawing();

            ClearBackground(WHITE);

            // Plot text
            DrawText(TextFormat("Time: %.2fs", (float)time), screenWidth * 0.05f,
                    screenHeight * 0.1f, 20, BLACK);
            DrawText(TextFormat("Potential Energy: %.2fJ",
                                (float)dip->getPotentialEnergy()),
                    screenWidth * 0.05f, screenHeight * 0.125f, 20, BLACK);
            DrawText(TextFormat("Kinetic Energy: %.2fJ", (float)dip->getKineticEnergy()),
                    screenWidth * 0.05f, screenHeight * 0.15f, 20, BLACK);
            DrawText(TextFormat("Total Energy: %.2fJ",
                                (float)(dip->getKineticEnergy() +
                                        dip->getPotentialEnergy())),
                    screenWidth * 0.05f, screenHeight * 0.175f, 20, BLACK);

            // Cart Track
            DrawLineEx((Vector2){screenWidth * 0.05, screenHeight * 0.55},
                    (Vector2){screenWidth * 0.7, screenHeight * 0.55}, 
                    5,
                    Fade(BLACK, 0.8f));
            // Cart (Darker with increasing mass)
            DrawRectangleV(cartPos, cartDims,
                        Fade(BLUE, initSliders[0].currValue / 20.0f + 0.5f));
            // Left wheel
            DrawCircleV(leftWheelPos, wheelRadius, Fade(RED, 0.95f));
            // Right wheel
            DrawCircleV(rightWheelPos, wheelRadius, Fade(RED, 0.95f));
            // Link 1
            DrawLineEx((Vector2){cartPos.x + cartDims.x / 2, cartPos.y + cartDims.y / 2},
                        massOnePos, 10, Fade(BLACK, 0.9f));
            // Link 2
            DrawLineEx(massOnePos, massTwoPos, 10, Fade(BLACK, 0.9f));
            // Mass 1 (Radius grows with mass)
            DrawCircleV(massOnePos, (sqrt(initSliders[1].currValue) + 1.0) * meter / 10, RED);
            // Mass 2 (Radius grows with mass)
            DrawCircleV(massTwoPos, (sqrt(initSliders[2].currValue) + 1.0) * meter / 10, RED);

            // Make GUI area on right quarter of window
            DrawLine(screenWidth * 0.75, 0, screenWidth * 0.75, screenWidth,
                    Fade(LIGHTGRAY, 0.6f));
            DrawRectangle(screenWidth * 0.75, 0, screenWidth * 0.25, screenHeight,
                        Fade(LIGHTGRAY, 0.3f));

            // Lock GUI sliders as they are not meant to accept input, but to show state
            GuiLock();
            // Add GUI sliders
            for (size_t i = 0; i < std::size(stateSliders); ++i) {
                GuiLabel((Rectangle){screenWidth * 0.775,
                                    screenHeight * (0.07f + i * 0.06f), 1400, 24},
                                    TextFormat(stateSliders[i].label, stateSliders[i].currValue));
                GuiSliderBar((Rectangle){screenWidth * 0.775,
                                        screenHeight * (0.09f + i * 0.06f),
                                        screenWidth * 0.2, screenHeight * 0.025},
                            NULL, NULL, &stateSliders[i].currValue, stateSliders[i].minValue,
                            stateSliders[i].maxValue);
            }

            // Control Mode
            GuiLabel(
                (Rectangle){screenWidth * 0.775, screenHeight * 0.01, 140, 24},
                "Control Mode");
            if (GuiDropdownBox(
                    (Rectangle){screenWidth * 0.775, screenHeight * 0.03,
                                screenWidth * 0.25 * 0.8, screenHeight * 0.025},
                    "Open Loop;Continuos LQR; Discrete LQR", &controlModeActive,
                    controlModeEdit))
                controlModeEdit = !controlModeEdit;

            GuiUnlock();

            // Pause/resume button
            if (pauseBtn)
                btnText = "Resume";
            else
                btnText = "Pause";
            if (GuiButton((Rectangle){screenWidth * 0.775, screenHeight * 0.43,
                                    screenWidth * 0.25 * 0.8, screenHeight * 0.025},
                        btnText)) {
                pauseBtn = !pauseBtn;
            }

            // Reset button
            if (GuiButton((Rectangle){screenWidth * 0.775, screenHeight * 0.47,
                                    screenWidth * 0.25 * 0.8, screenHeight * 0.025},
                        "Reset")) {
                resetBtn = !resetBtn;
            }

            EndDrawing();
        }
        
        // Reset button was pressed, reset buttons and time for next sim
        if (resetBtn) {
            resetBtn = false;
            startBtn = false;
            pauseBtn = false;
            time = 0;
            delete dip;
            dip = nullptr;
            delete lqr;
            lqr = nullptr;
        }
    }

    // Delete dip;
    delete dip;

    CloseWindow();

    return 0;
}