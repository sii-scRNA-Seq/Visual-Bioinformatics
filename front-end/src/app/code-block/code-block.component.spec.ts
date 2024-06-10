import { By } from '@angular/platform-browser';
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { FormsModule } from '@angular/forms';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatCardModule } from '@angular/material/card';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { BlockService } from '../block.service';
import { CodeBlockComponent } from './code-block.component';
import { MockBlockService } from '../mock-block.service';
import { OutputService } from '../output.service';
import { MockOutputService } from '../mock-output.service';

describe('CodeBlockComponent', () => {
  let component: CodeBlockComponent;
  let fixture: ComponentFixture<CodeBlockComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [CodeBlockComponent],
      imports: [
        HttpClientTestingModule,
        MatCardModule,
        MatSnackBarModule,
        FormsModule
      ],
      providers: [
        { provide: BlockService, useClass: MockBlockService },
        { provide: OutputService, useClass: MockOutputService },
      ],
    });
    fixture = TestBed.createComponent(CodeBlockComponent);
    component = fixture.componentInstance;
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('Removing Blocks', () => {
    it('blockService.removeBlock should be called with loaddata when remove button is clicked', () => {
      component.block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('loaddata');
    });

    it('blockService.removeBlock should be called with basicfiltering when remove button is clicked', () => {
      component.block = {
        blockId: 'basicfiltering',
        title: 'Basic Filtering',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('basicfiltering');
    });

    it('blockService.removeBlock should be called with qcplots when remove button is clicked', () => {
      component.block = {
        blockId: 'qcplots',
        title: 'Quality Control Plots',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('qcplots');
    });

    it('blockService.removeBlock should be called with qcfiltering when remove button is clicked', () => {
      component.block = {
        blockId: 'qcfiltering',
        title: 'Quality Control Filtering',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('qcfiltering');
    });

    it('blockService.removeBlock should be called with variablegenes when remove button is clicked', () => {
      component.block = {
        blockId: 'variablegenes',
        title: 'Identify Highly Variable Genes',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('variablegenes');
    });

    it('blockService.removeBlock should be called with pca when remove button is clicked', () => {
      component.block = {
        blockId: 'pca',
        title: 'Principal Component Analysis',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('pca');
    });

    it('blockService.removeBlock should be called with runumap when remove button is clicked', () => {
      component.block = {
        blockId: 'runumap',
        title: 'Run UMAP',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('runumap');
    });

    it('should be available when blocks are not being executed', () => {
      component.block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      component.executingBlocks = false;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(false);
    });

    it('should become disabled while blocks are being executed', () => {
      component.block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      component.executingBlocks = false;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(false);
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(true);
    });

    it('should become available once blocks have stopped being executed', () => {
      component.block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(true);
      component.executingBlocks = false;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(false);
    });
  });

  describe('Parameter Inputs', () => {
    it('should be available when blocks are not being executed', () => {
      component.block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {type: 'InputParameter', key: 'test_param', text: 'Test Parameter', value: 0},
        ],
      };
      fixture.detectChanges();
      component.executingBlocks = false;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('input')).nativeElement.disabled).toBe(false);
    });

    it('should become disabled while blocks are being executed', () => {
      component.block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {type: 'InputParameter', key: 'test_param', text: 'Test Parameter', value: 0},
        ],
      };
      fixture.detectChanges();
      component.executingBlocks = false;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('input')).nativeElement.disabled).toBe(false);
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('input')).nativeElement.disabled).toBe(true);
    });

    it('should become available once blocks have stopped being executed', () => {
      component.block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {type: 'InputParameter', key: 'test_param', text: 'Test Parameter', value: 0},
        ],
      };
      fixture.detectChanges();
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('input')).nativeElement.disabled).toBe(true);
      component.executingBlocks = false;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('input')).nativeElement.disabled).toBe(false);
    });
  });
});
