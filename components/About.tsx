
import React from 'react';

const About: React.FC = () => {
  return (
    <section id="about" className="py-20 md:py-32">
      <h2 className="text-3xl md:text-4xl font-bold text-center text-stone-100 mb-12">About Me</h2>
      <div className="flex flex-col md:flex-row items-center gap-12">
        <div className="md:w-1/3 flex-shrink-0">
          <img 
            src="https://picsum.photos/seed/profile/400/400" 
            alt="Profile portrait" 
            className="rounded-full shadow-lg w-64 h-64 md:w-80 md:h-80 mx-auto object-cover border-4 border-lime-500/50"
          />
        </div>
        <div className="md:w-2/3 text-lg text-stone-400 space-y-4 text-center md:text-left">
          <p>
            Hello! I'm a passionate Senior Frontend Engineer with a deep expertise in creating modern, responsive, and user-friendly web applications. With over a decade of experience in the field, I specialize in the React ecosystem, TypeScript, and state-of-the-art UI/UX design principles.
          </p>
          <p>
            My journey in web development began with a fascination for how code could be transformed into interactive and engaging experiences. This passion drives me to constantly learn and adapt to new technologies, ensuring that the products I build are not only functional but also delightful to use.
          </p>
          <p>
            When I'm not coding, you can find me exploring new hiking trails, contributing to open-source projects, or diving into a good book on design theory.
          </p>
        </div>
      </div>
    </section>
  );
};

export default About;