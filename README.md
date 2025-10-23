# React Portfolio for GitHub Pages (Vite Edition)

This is a modern, single-page portfolio application built with React, TypeScript, and Tailwind CSS, and powered by the Vite build tool. It's structured for a professional development workflow and easy deployment to static hosts like GitHub Pages.

## Why the Old Method Didn't Work (The "Green Screen" Problem)

You were seeing a blank green screen because web browsers cannot directly run TypeScript (`.tsx`) or JSX code. They only understand plain JavaScript.

Your project was missing a **build step**, which is a crucial process that does two things:
1.  **Transpiles**: Converts your TypeScript/JSX code into browser-readable JavaScript.
2.  **Bundles**: Packages all your code and dependencies into a few optimized files for efficient loading.

Opening the `index.html` file locally (`file://...`) also fails because modern JavaScript modules have security restrictions that prevent them from working without a web server.

This updated project structure uses **Vite** to provide both a local development server and a build process, solving all these issues.

---

## How to Deploy to GitHub Pages (The Correct Way)

Follow these steps to get your portfolio online.

### Prerequisites

You need to have [Node.js](https://nodejs.org/) installed on your computer. This will also install `npm`.

### Step 1: Set Up the Project Locally

1.  **Download the Code**: Download all the files for this project to a folder on your computer.
2.  **Open a Terminal**: Open your terminal or command prompt and navigate into the project folder.
3.  **Install Dependencies**: Run the following command to install all the necessary packages (React, Vite, etc.):
    ```bash
    npm install
    ```

### Step 2: Run the Local Development Server

To see your website and test changes locally, run:
```bash
npm run dev
```
This will start a local server, and you can view your portfolio at the URL it provides (usually `http://localhost:5173`). The site will automatically reload as you save changes to the code.

### Step 3: Create a GitHub Repository

1.  Go to [GitHub](https://github.com) and create a new public repository named `your-username.github.io` (replace `your-username` with your GitHub username).
2.  Follow GitHub's instructions to push your local project folder to this new repository.

### Step 4: Build the Project for Deployment

Before you can deploy, you need to create the optimized production-ready files. In your terminal, run:
```bash
npm run build
```
This command will create a new `dist` folder in your project. This `dist` folder contains the final, static HTML, CSS, and JavaScript files that will be hosted.

### Step 5: Deploy to GitHub Pages

We will use the `gh-pages` package to make deployment simple.

1.  **Install `gh-pages`**: In your terminal, run this command:
    ```bash
    npm install gh-pages --save-dev
    ```
2.  **Add a Deploy Script**: Open your `package.json` file and add a `"deploy"` script inside the `"scripts"` section:
    ```json
    "scripts": {
      "dev": "vite",
      "build": "tsc && vite build",
      "deploy": "gh-pages -d dist", // Add this line
      "preview": "vite preview"
    },
    ```
3.  **Add Homepage URL**: Still in `package.json`, add a new `homepage` field at the top level, specifying your GitHub Pages URL:
    ```json
    {
      "name": "react-portfolio",
      "private": true,
      "version": "0.0.0",
      "type": "module",
      "homepage": "https://your-username.github.io", // Add this line
      // ... rest of the file
    }
    ```
    **Important**: Replace `your-username` with your actual GitHub username.

4.  **Run the Deploy Script**: Now, run the deploy command in your terminal:
    ```bash
    npm run deploy
    ```
This will automatically push the contents of your `dist` folder to a special `gh-pages` branch on your repository.

5.  **Configure GitHub Pages Source**:
    *   In your GitHub repository, go to **Settings > Pages**.
    *   Under "Build and deployment", set the **Source** to **Deploy from a branch**.
    *   Set the **Branch** to `gh-pages` and the folder to `/ (root)`. Click **Save**.

### Step 6: You're Live!

Your portfolio should be live at `https://your-username.github.io` within a few minutes.
