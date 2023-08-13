using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FluidSim2 : MonoBehaviour
{

	static int N = 64;
    static int size = (N + 2) * (N + 2);
    float[] u, v, u_prev, v_prev, dens, dens0;
    [SerializeField] float dt = 0.01f;
    [SerializeField] float diff = 0f;
    [SerializeField] float visc = 0f;

    private Texture2D texture;
    SpriteRenderer m_SpriteRenderer;
    Color[] colors;

    // Start is called before the first frame update
    void Start()
    {
        u = new float[size];                // x
        v = new float[size];                // y
        u_prev = new float[size];
        v_prev = new float[size];
        dens = new float[size];                //density
        dens0 = new float[size];

        texture = new Texture2D(N+2, N+2, TextureFormat.RGBAFloat, -1, false);
        m_SpriteRenderer = GetComponent<SpriteRenderer>();
        colors = new Color[size];

        MaterialPropertyBlock block = new MaterialPropertyBlock();
        block.SetTexture("_MainTex", texture);
        m_SpriteRenderer.SetPropertyBlock(block);

        add_source();
    }

    // Update is called once per frame
    void Update()
    {
        //add_source2();
        velStep(N, u, v, u_prev, v_prev, visc, dt);
        densStep(N, dens, dens0, u, v, diff, dt);


        updateColors();
    }

    private void updateColors()
    {
        for (int c = 0; c < size; c++)
        {
            colors[c].r = dens[c];
            colors[c].g = 0;
                colors[c].b = 0;
            colors[c].a = 1;
        }

        texture.SetPixels(colors);
        texture.Apply();
    }

    private int Cell(int i, int j)
    {
        return (i + ((N + 2) * j));
    }

    void add_source()
    {
        for (int c = 500; c < size-4; c++)
        {
            dens[c] = Random.Range(0f, 0.5f);
        }

    }


    /**
     * Diffuse a field at rate diff as per Stam(2003). Iterative solver for stable diffusion. 
     * 
    **/
    void diffuse(int n, int b, float[] current_field, float[] prev_field, float diff, float dt)
    {
        int x, y, k;
        float a = dt * diff * n * n;

        for (k = 0; k < 20; k++)
        {
            for (x = 1; x <= n; x++)
            {
                for (y = 1; y <= n; y++)
                {
                    current_field[Cell(x, y)] = (prev_field[Cell(x, y)] +  a * (current_field[Cell(x - 1, y)] + current_field[Cell(x + y, y)] + current_field[Cell(x, y - 1)] + current_field[Cell(x, y + 1)])) / (1 + (4 * a));
                }
            }
            setBnd(n, b, current_field);
        }
    }


    void advect(int n, int b, float[] current_field, float[] prev_field, float[] u, float[] v, float dt)
    {
        int i, j, i0, j0, i1, j1;
        float x, y, s0, t0, s1, t1, dt0;

        dt0 = dt * n;
        for (i = 1; i <= n; i++)
        {
            for (j = 1; j <= n; j++)
            {
                x = i - dt0 * u[Cell(i, j)];
                y = j - dt0 * v[Cell(i, j)];

                if (x < 0.5f)
                {
                    x = 0.5f;
                }
                if (x > (n + 0.5f))
                {
                    x = n + 0.5f;
                }
                i0 = (int)x;
                i1 = i0 + 1;

                if (y < 0.5f)
                {
                    y = 0.5f;
                }
                if (y > (n + 0.5f))
                {
                    y = n + 0.5f;
                }
                j0 = (int)y;
                j1 = j0 + 1;

                s1 = x - i0;
                s0 = 1 - s1;

                t1 = y - j0;
                t0 = 1 - t1;

                current_field[Cell(i, j)] = s0 * (t0 * prev_field[Cell(i0, j0)] + t1 * prev_field[Cell(i0, j1)]) +
                  s1 * (t0 * prev_field[Cell(i1, j0)] + t1 * prev_field[Cell(i1, j1)]);
            }
        }
        setBnd(n, b, current_field);
    }


    void setBnd(int n, int b, float[] x)
    {
        int i;

        for (i = 0; i <= n; i++)
        {
            x[Cell(0, i)] = (b == 1 ? -x[Cell(1, i)] : x[Cell(1, i)]);
            x[Cell(n + 1, i)] = (b == 1 ? -x[Cell(n, i)] : x[Cell(n, i)]);
            x[Cell(i, 0)] = (b == 2 ? -x[Cell(i, 1)] : x[Cell(i, 1)]);
            x[Cell(i, n + 1)] = (b == 2 ? -x[Cell(i, n)] : x[Cell(i, n)]);
        }
        x[Cell(0, 0)] = 0.5f * (x[Cell(1, 0)] + x[Cell(0, 1)]);
        x[Cell(0, n + 1)] = 0.5f * (x[Cell(1, n + 1)] + x[Cell(0, n)]);
        x[Cell(n + 1, 0)] = 0.5f * (x[Cell(n, 0)] + x[Cell(n + 1, 1)]);
        x[Cell(n + 1, n + 1)] = 0.5f * (x[Cell(n, n + 1)] + x[Cell(n + 1, n)]);
    }


    void densStep(int n, float[] x, float[] x0, float[] u, float[] v, float diff, float dt)
    {
        float[] temp = x0;
        x0 = x;
        x = temp;
        diffuse(n, 0, x, x0, diff, dt);
        temp = x0;
        x0 = x;
        x = temp;
        advect(n, 0, x, x0, u, v, dt);
    }


    void velStep(int n,
    float[] u,
    float[] v,
    float[] u0,
    float[] v0,
    float visc,
    float dt)
    {
        //    swap(u0, u);
        float[] temp = u0;
        u0 = u;
        u = temp;

        diffuse(n, 1, u, u0, visc, dt);

        //    swap(v0, v);
        temp = v0;
        v0 = v;
        v = temp;

        diffuse(n, 2, v, v0, visc, dt);
        project(n, u, v, u0, v0);

        //    swap(u0, u);
        temp = u0;
        u0 = u;
        u = temp;

        //    swap(v0, v);
        temp = v0;
        v0 = v;
        v = temp;

        advect(n, 1, u, u0, u0, v0, dt);
        advect(n, 2, v, v0, u0, v0, dt);
        project(n, u, v, u0, v0);
       
    }

    void project(int n,
        float[] u,
        float[] v,
        float[] p,
        float[] div)
    {
        int i, j, k;
        float h = 1.0f / n;
        for (i = 1; i <= n; i++)
        {
            for (j = 1; j <= n; j++)
            {
                div[Cell(i, j)] = -0.5f * h * (u[Cell(i + 1, j)] - u[Cell(i - 1, j)] +
                    v[Cell(i, j + 1)] - v[Cell(i, j - 1)]);
                p[Cell(i, j)] = 0;
            }
        }
        setBnd(n, 0, div);
        setBnd(n, 0, p);

        for (k = 0; k < 20; k++)
        {
            for (i = 1; i <= n; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    p[Cell(i, j)] = (div[Cell(i, j)] + p[Cell(i - 1, j)] + p[Cell(i + 1, j)] +
                        p[Cell(i, j - 1)] + p[Cell(i, j + 1)]) / 4;
                }
            }
            setBnd(n, 0, p);
        }

        for (i = 1; i <= n; i++)
        {
            for (j = 1; j <= n; j++)
            {
                u[Cell(i, j)] -= 0.5f * (p[Cell(i + 1, j)] - p[Cell(i - 1, j)]) / h;
                v[Cell(i, j)] -= 0.5f * (p[Cell(i, j + 1)] - p[Cell(i, j - 1)]) / h;
            }
        }
        setBnd(n, 1, u);
        setBnd(n, 2, v);
    }

}
