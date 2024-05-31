#include "raylib.h"

int main(void)
{
    InitWindow(1600, 900, "DIP-Control");

    while (!WindowShouldClose())
    {
        BeginDrawing();
            ClearBackground(RAYWHITE);
            DrawText("DIP-Control", 190, 200, 20, LIGHTGRAY);
        EndDrawing();
    }

    CloseWindow();

    return 0;
}