import { ComponentFixture, TestBed } from '@angular/core/testing';
import { By } from '@angular/platform-browser';
import { MatCardModule } from '@angular/material/card';
import { from } from 'rxjs';
import { CodeBlockComponent } from './code-block.component';
import { BlockService } from '../block.service';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatSnackBarModule } from '@angular/material/snack-bar';

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
      ],
      providers: [BlockService],
    });
    fixture = TestBed.createComponent(CodeBlockComponent);
    component = fixture.componentInstance;
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('RemoveBlock', () => {
    it ('should be called with loaddata when Load Data block is removed', () => {
      component.block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
        onRun: () => from(''),
      };
      fixture.detectChanges();

      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('loaddata');
    });

    it ('should be called with basicfiltering when Basic Filtering block is removed', () => {
      component.block = {
        blockId: 'basicfiltering',
        title: 'Basic Filtering',
        possibleChildBlocks: [],
        parameters: [],
        onRun: () => from(''),
      };
      fixture.detectChanges();

      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('basicfiltering');
    });
  });
});
